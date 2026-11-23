"""Triage main script: one pass over the undetermined reads, several reports."""

import argparse
import csv
import logging
import pathlib
from collections import Counter
from multiprocessing import Pool

from biomate.headers import count_header_indexes
from biomate.index.index import group_indexes
from biomate.setup import setup_logging
from biomate.strainer import strainer
from biomate.triage import barcodes
from biomate.waif import waif

# Most frequent indexes grouped per file, and groups reported per file
TOP_INDEXES = 1000
TOP_GROUPS = 20


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    description = (
        "Triage the undetermined reads of a flowcell: strainer, index "
        "discovery and waif from a single pass over the FASTQ headers."
    )
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=description,
    )
    strainer.add_arguments(parser)
    parser.set_defaults(parse=strainer.validate_args, run=main)

    return parser


def _scan(path: pathlib.Path) -> tuple[int, list[dict], list[tuple[str, int]]]:
    """Count one file; return only what the reports need (workers pickle it)."""
    # Serial count: the threaded inflate of headers.py measured 2.5x slower
    # on a 4.7 GB NextSeq undetermined file (234 s against 90 s)
    counts = count_header_indexes(path)
    return (
        sum(counts.values()),
        strainer.top_undetermined(path, counts),
        counts.most_common(TOP_INDEXES),
    )


def write_top_indexes(
    scans: list[tuple[pathlib.Path, int, list[tuple[str, int]]]],
    input_path: pathlib.Path,
    output_file: pathlib.Path,
) -> None:
    """Write the largest index groups (1 substitution) of every file."""
    with open(output_file, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["File", "Index", "Count", "Fraction", "Variants"])
        for path, total, top in scans:
            groups = group_indexes(Counter(dict(top)), 1)
            sizes = sorted(
                (
                    (sum(c for _, _, c in members), seed)
                    for seed, members in groups.items()
                ),
                key=lambda item: (-item[0], item[1]),
            )
            for count, seed in sizes[:TOP_GROUPS]:
                writer.writerow(
                    [
                        path.relative_to(input_path),
                        seed,
                        count,
                        f"{count / total:.4f}",
                        len(groups[seed]),
                    ]
                )


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running triage module...")
    basename = args.input_path.resolve().name

    sample_sheet = args.sample_sheet or args.input_path / "SampleSheet.csv"
    sample_sheet_indexes = strainer.extract_indexes_from_sample_sheet(sample_sheet)

    files = strainer.find_undetermined(args.input_path)
    logging.info(f"Counting the header indexes of {len(files)} undetermined files...")
    if args.threads > 1 and len(files) > 1:
        with Pool(min(args.threads, len(files))) as pool:
            scans = pool.map(_scan, files)
    else:
        scans = [_scan(f) for f in files]

    logging.info("Strainer: searching for unexpected indexes...")
    strainer.report(
        sample_sheet_indexes,
        [row for _, rows, _ in scans for row in rows],
        args.output_path,
        basename,
    )

    logging.info("Barcode findings: explaining undetermined indexes...")
    barcodes.report(
        sample_sheet,
        args.input_path,
        [row for _, rows, _ in scans for row in rows],
        args.output_path,
        basename,
    )

    top_file = args.output_path / f"{basename}_top_undetermined.csv"
    write_top_indexes(
        [(f, total, top) for f, (total, _, top) in zip(files, scans) if total],
        args.input_path,
        top_file,
    )
    logging.info(f"Top undetermined indexes saved to {top_file}.")

    # Waif needs re-run datasets (Stats.json); it reads no FASTQ file here
    try:
        table = waif.load_stats_frames(args.input_path)
    except ValueError as e:
        logging.info(f"Waif skipped: {e}")
        return
    lanes = waif.find_shared_lanes(table)
    if not lanes:
        logging.info("Waif skipped: no lane appears in more than one dataset.")
        return
    waif_file = args.output_path / f"{basename}_always_undetermined.csv"
    waif.shared_undetermined(table, lanes).write_csv(waif_file)
    logging.info(f"Always-undetermined barcodes saved to {waif_file}.")
