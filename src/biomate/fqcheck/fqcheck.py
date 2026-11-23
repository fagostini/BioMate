"""FqCheck main script: validate FASTQ files and their mates before delivery."""

import argparse
import csv
import hashlib
import logging
import pathlib
import re
import sys
import zlib
from multiprocessing import Pool

import dnaio

from biomate.setup import setup_logging

# Read token of an Illumina file name, e.g. Sample_S1_L001_R1_001.fastq.gz
READ_TOKEN = re.compile(r"_(R1|R2|I1|I2)_")
ILLUMINA_NAME = re.compile(r"_S\d+(_L\d{3})?_[RI][12]_001\.fastq(\.gz)?$")
# Mates are read in this order, so R1 always leads the pairing check
READ_ORDER = ("R1", "R2", "I1", "I2")


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    description = (
        "Check that FASTQ files are complete and that their mates (R1, R2, "
        "I1, I2) hold the same reads, optionally writing md5 checksums."
    )
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=description,
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        required=True,
        help="Directory searched recursively for *.fastq.gz and *.fastq files.",
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        default=pathlib.Path("./"),
        help="Directory where the report is written (default: current directory).",
    )
    parser.add_argument(
        "--md5",
        action="store_true",
        help="Also write an md5sum-compatible checksum file.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of file groups checked in parallel.",
    )
    parser.set_defaults(parse=validate_args, run=main)
    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} is not a valid directory."
        )
    if args.threads < 1:
        raise argparse.ArgumentTypeError("--threads must be >= 1.")
    if not args.output_path.is_dir():
        logging.info(
            f"Output directory {args.output_path} does not exist. Creating it..."
        )
        args.output_path.mkdir(parents=True, exist_ok=True)
    return args


def find_fastqs(input_path: pathlib.Path) -> list[pathlib.Path]:
    """FASTQ files under ``input_path``; a link and its target count once."""
    files: dict[pathlib.Path, pathlib.Path] = {}
    for pattern in ("**/*.fastq.gz", "**/*.fastq"):
        for f in sorted(input_path.glob(pattern)):
            if not f.is_file():
                logging.warning(f"Skipping {f}: broken link")
                continue
            files.setdefault(f.resolve(), f)
    return sorted(files.values())


def group_mates(files: list[pathlib.Path]) -> dict[str, dict[str, pathlib.Path]]:
    """Group files by name with the read token masked: {group: {read: file}}.

    A file without a read token is a group of its own, read as R1.
    """
    groups: dict[str, dict[str, pathlib.Path]] = {}
    for f in files:
        match = READ_TOKEN.search(f.name)
        read = match.group(1) if match else "R1"
        key = str(f.with_name(READ_TOKEN.sub("_*_", f.name, count=1)))
        groups.setdefault(key, {})[read] = f
    return groups


def check_group(mates: dict[str, pathlib.Path], md5: bool) -> dict:
    """Read every mate once in lockstep; return the group's row and checksums."""
    files = [mates[r] for r in READ_ORDER if r in mates]
    warnings = [
        f"{f.name}: not an Illumina file name"
        for f in files
        if not ILLUMINA_NAME.search(f.name)
    ]
    result = {"files": files, "records": 0, "status": "ok", "message": ""}
    if "R1" not in mates:
        result.update(status="error", message="no R1 file for this group")
    else:
        try:
            with dnaio.open(*files) as reader:
                result["records"] = sum(1 for _ in reader)
        except (dnaio.FileFormatError, EOFError, OSError, zlib.error) as e:
            result.update(status="error", message=f"{type(e).__name__}: {e}")
        else:
            if result["records"] == 0:
                warnings.append("no records")
    if result["status"] == "ok" and warnings:
        result.update(status="warning", message="; ".join(warnings))
    result["md5"] = {}
    if md5:
        for f in files:
            with open(f, "rb") as fh:
                result["md5"][f] = hashlib.file_digest(fh, "md5").hexdigest()
    return result


def _check(item: tuple[dict[str, pathlib.Path], bool]) -> dict:
    return check_group(*item)


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running fqcheck module...")
    basename = args.input_path.resolve().name

    groups = group_mates(find_fastqs(args.input_path))
    if not groups:
        raise ValueError(f"No FASTQ files found in {args.input_path}.")
    logging.info(f"Checking {len(groups)} FASTQ file groups...")
    items = [(mates, args.md5) for mates in groups.values()]
    if args.threads > 1 and len(items) > 1:
        with Pool(min(args.threads, len(items))) as pool:
            results = pool.map(_check, items)
    else:
        results = [_check(item) for item in items]

    report = args.output_path / f"{basename}_fqcheck.csv"
    with open(report, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Group", "Files", "Records", "Status", "Message"])
        for key, r in zip(groups, results):
            writer.writerow(
                [
                    pathlib.Path(key).relative_to(args.input_path),
                    ";".join(f.name for f in r["files"]),
                    r["records"],
                    r["status"],
                    r["message"],
                ]
            )
    logging.info(f"Report saved to {report}.")

    if args.md5:
        md5_file = args.output_path / f"{basename}.md5"
        with open(md5_file, "w") as handle:
            for r in results:
                for f, digest in r["md5"].items():
                    handle.write(f"{digest}  {f.relative_to(args.input_path)}\n")
        logging.info(f"Checksums saved to {md5_file}.")

    errors = sum(r["status"] == "error" for r in results)
    warnings = sum(r["status"] == "warning" for r in results)
    logging.info(f"{len(results)} groups: {errors} error(s), {warnings} warning(s).")
    if errors:
        sys.exit(1)
