"""Waif main script.

Report barcodes that remain undetermined in every demultiplexing run of a
lane.

For each lane that appears in two or more datasets (bcl-convert runs stored in
``<input_path>/<Dataset>/Stats/Stats.json`` directories, or in
``<input_path>/<Dataset>/Reports/Top_Unknown_Barcodes.csv`` with
``--bclconvert-reports``), the barcodes that
are undetermined in *all* of those datasets are reported with an estimated
read count: the mean of the datasets' counts, rounded up to a multiple of 10.
"""

import argparse
import logging
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor

import pathlib

import polars

from biomate.headers import _count_for_file
from biomate.samplesheet import read_section
from biomate.setup import setup_logging

STATS_FILE_NAME = "Stats.json"
TOP_UNKNOWN_FILE_NAME = "Top_Unknown_Barcodes.csv"
EXCLUDED_DIR = "Demultiplexing"
INDEX_SEP = "+"
COUNT_ROUND_TO = 10
COUNT_COLUMN_PREFIX = "Count_"


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description=(
            "Report barcodes that remain undetermined in every demultiplexing "
            "run of a lane."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=(
            "Report barcodes that remain undetermined in every demultiplexing "
            "run of a lane."
        ),
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        required=True,
        help=(
            "Path to the input directory: the run output with Stats.json files "
            "or, with --exact, the undetermined FASTQ files."
        ),
    )
    parser.add_argument(
        "--output-file",
        type=pathlib.Path,
        required=False,
        help=(
            "File to write the results to. If not specified, the results are "
            "written to stdout as CSV."
        ),
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help=(
            "Enable exact counting of undetermined indexes by reading the "
            "undetermined FASTQ files. Needs a Sample Sheet."
        ),
    )
    parser.add_argument(
        "--bclconvert-reports",
        action="store_true",
        help=(
            "Read the bcl-convert Reports/Top_Unknown_Barcodes.csv files "
            "instead of Stats.json."
        ),
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        required=False,
        help=(
            "Path to the Sample Sheet file (used with --exact). If not "
            "provided, the module looks for SampleSheet.csv in the input "
            "directory."
        ),
    )
    parser.add_argument(
        "--lanes",
        type=int,
        nargs="+",
        default=None,
        metavar="LANE",
        help=(
            "Only process these lanes (e.g. --lanes 1 3). Lanes that would "
            "not be selected anyway (not shared between datasets, or not "
            "multi-project in the Sample Sheet with --exact) are skipped "
            "with a warning. Default: all selected lanes."
        ),
    )
    parser.add_argument(
        "--block-cache",
        type=pathlib.Path,
        default=None,
        metavar="DIR",
        help=(
            "Directory for a small gzip member-offset index of the "
            "undetermined files; enables forked per-file workers that "
            "decompress and count in parallel (effective when there are "
            "fewer files than cores). Each file is scanned once to locate "
            "member boundaries (a few minutes per large file); the index "
            "file is tiny and the source is never copied or re-compressed. "
            "Only used with --exact; Linux only."
        ),
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} does not exist or is not a directory."
        )
    if args.output_file is not None and not args.output_file.parent.is_dir():
        logging.info(
            f"Output directory {args.output_file.parent} does not exist. Creating it..."
        )
        args.output_file.parent.mkdir(parents=True, exist_ok=True)
    if args.exact and args.bclconvert_reports:
        raise argparse.ArgumentTypeError(
            "--exact and --bclconvert-reports cannot be used together."
        )
    if args.exact:
        if args.sample_sheet is not None and not args.sample_sheet.is_file():
            raise argparse.ArgumentTypeError(
                f"Sample Sheet {args.sample_sheet} does not exist."
            )
        sample_sheet = args.sample_sheet or args.input_path / "SampleSheet.csv"
        if not sample_sheet.is_file():
            raise argparse.ArgumentTypeError(
                f"Sample Sheet not found in {args.input_path} "
                "(required for --exact); provide one with --sample-sheet."
            )
    if args.lanes is not None and any(lane < 1 for lane in args.lanes):
        raise argparse.ArgumentTypeError(
            f"Lane numbers must be positive integers, got {args.lanes}."
        )
    if args.block_cache is not None:
        if sys.platform != "linux":
            raise argparse.ArgumentTypeError(
                "--block-cache is only supported on Linux (it relies on os.fork)."
            )
        args.block_cache.mkdir(parents=True, exist_ok=True)
    return args


def get_shared_lanes_from_samplesheet(
    samplesheet: pathlib.Path,
) -> list[int]:
    """Return the lanes carrying samples from two or more Sample_Projects.

    Rows without any index (``index``/``index2`` both empty) and rows with an
    empty project are ignored.
    """
    projects_by_lane: defaultdict[str, set[str]] = defaultdict(set)
    data = read_section(samplesheet, infer_schema=False, bare_header=("Lane", "FCID"))
    columns = ["Lane", "index", "index2", "Sample_Project"]
    data = data.with_columns(
        polars.lit(None, dtype=polars.String).alias(col)
        for col in columns
        if col not in data.columns
    )
    for lane, index1, index2, project in data.select(columns).iter_rows():
        if (index1 or index2) and project:
            projects_by_lane[lane].add(project)

    return sorted(
        int(lane) for lane, projects in projects_by_lane.items() if len(projects) > 1
    )


def _count_file(
    task: tuple[int, pathlib.Path, int, pathlib.Path | None],
) -> tuple[int, pathlib.Path, dict[str, int]]:
    """Process-pool worker: count one FASTQ file (module-level so it pickles)."""
    lane, f, budget, cache_dir = task
    return lane, f, dict(_count_for_file(f, budget, cache_dir))


def _barcode_frame(counts: Counter, lane: int, f: pathlib.Path) -> polars.DataFrame:
    """Build the long barcode table of one counted file."""
    return (
        polars.DataFrame(
            {"Barcode": list(counts), "Count": list(counts.values())},
            schema={"Barcode": polars.String, "Count": polars.Int64},
        )
        .with_columns(
            polars.lit(lane, dtype=polars.Int64).alias("Lane"),
            polars.lit(f.parent.name).alias("Dataset"),
            polars.col("Barcode").str.len_chars().alias("Barcode_len"),
        )
        .select(["Lane", "Barcode", "Count", "Dataset", "Barcode_len"])
        .with_columns(
            polars.col("Barcode")
            .str.split_exact(INDEX_SEP, 1)
            .struct.rename_fields(["index1", "index2"])
            .alias("indexes")
        )
        .unnest("indexes")
        # A single-index barcode (no separator) yields a null index2; normalize
        # it to "" so a dataset is either uniformly single (index2 == "") or
        # uniformly dual (index2 non-empty).
        .with_columns(polars.col("index2").fill_null(""))
    )


def fastqs_to_stats_frames(
    input_dir: pathlib.Path,
    selected_lanes: list[int],
    block_cache: pathlib.Path | None = None,
) -> polars.DataFrame:
    """Convert the selected undetermined FASTQ files to a long barcode table.

    Expects the FASTQ files directly under ``<Dataset>/`` (the dataset name is
    ``file.parent.name``); only R1 files are used, as the index information is
    in the R1 headers. Files below any ``Demultiplexing/`` directory are
    excluded (that is the primary bcl-convert demux output, not a re-run). If
    bcl-convert split a dataset into several 1M-read parts (``_001``, ``_002``,
    ...), their counts are summed, so each (Lane, Dataset, Barcode) row is
    unique. Per-file counting runs in a process pool (one worker per file, up
    to the core count); the result order matches the file order, so the
    output is deterministic. Each file gets a core budget of
    ``cpu_count // len(files)``: with enough slack and a ``block_cache``
    directory it is counted by forked workers over a gzip member-offset index
    (a single-member file falls back to the serial pass); with slack but no
    ``block_cache`` it uses the threaded inflate+count pipeline; otherwise a
    single fused serial pass.

    Raises:
        ValueError: If no undetermined FASTQ files are found for the
            selected lanes.
    """
    lanes = set(selected_lanes)
    files = [
        (lane, f)
        for lane in sorted(lanes)
        for f in sorted(input_dir.glob(f"**/Undetermined_*_L{lane:03d}_R1_*.fastq.gz"))
        if EXCLUDED_DIR not in f.parts
    ]
    if not files:
        raise ValueError(
            f"No undetermined FASTQ files found for lanes {sorted(lanes)} "
            f"under {input_dir}."
        )
    for lane in sorted(lanes - {lane for lane, _ in files}):
        logging.warning(f"No undetermined FASTQ files for lane {lane}; skipping")

    cpu = os.cpu_count() or 1
    budget = max(1, cpu // len(files))
    tasks = [(lane, f, budget, block_cache) for lane, f in files]
    workers = min(len(files), cpu)
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_count_file, tasks, chunksize=1))
    else:
        results = [_count_file(task) for task in tasks]
    frames = [_barcode_frame(counts, lane, f) for lane, f, counts in results]
    return (
        polars.concat(frames)
        .group_by(["Lane", "Dataset", "Barcode"])
        .agg(
            polars.sum("Count").alias("Count"),
            polars.col("Barcode_len").first().alias("Barcode_len"),
            polars.col("index1").first().alias("index1"),
            polars.col("index2").first().alias("index2"),
        )
    )


def load_stats_frames(input_dir: pathlib.Path) -> polars.DataFrame:
    """Load every Stats.json under ``input_dir`` into a long barcode table.

    Expected layout: ``input_dir/<Dataset>/Stats/Stats.json``; the dataset
    name is the directory that contains the ``Stats`` folder
    (``file.parents[1]``). Files below any ``Demultiplexing/`` directory are
    excluded (that is the primary bcl-convert demux output, not a re-run).

    Raises:
        ValueError: If no Stats.json files are found, or if none of them
            carries an ``UnknownBarcodes`` field.
    """
    stats_files = [
        f
        for f in input_dir.glob(f"**/{STATS_FILE_NAME}")
        if EXCLUDED_DIR not in f.parts
    ]
    if not stats_files:
        raise ValueError(f"No {STATS_FILE_NAME} files found under {input_dir}.")

    frames = []
    for f in stats_files:
        try:
            frame = polars.read_json(f).select("UnknownBarcodes")
        except polars.exceptions.ColumnNotFoundError:
            logging.warning(f"{f} has no 'UnknownBarcodes' field; skipping")
            continue
        frame = (
            frame.explode("UnknownBarcodes", empty_as_null=True)
            .unnest("UnknownBarcodes")
            .unnest("Barcodes")
            .unpivot(index="Lane", variable_name="Barcode", value_name="Count")
            .drop_nulls()
            .with_columns(
                polars.lit(f.parents[1].name).alias("Dataset"),
                polars.col("Barcode").str.len_chars().alias("Barcode_len"),
            )
            .with_columns(
                polars.col("Barcode")
                .str.split_exact(INDEX_SEP, 1)
                .struct.rename_fields(["index1", "index2"])
                .alias("indexes")
            )
            .unnest("indexes")
            # A single-index barcode (no separator) yields a null index2;
            # normalize it to "" (see _barcode_frame).
            .with_columns(polars.col("index2").fill_null(""))
        )
        frames.append(frame)

    if not frames:
        raise ValueError(
            f"No usable {STATS_FILE_NAME} files found (all lacked 'UnknownBarcodes')."
        )
    return polars.concat(frames)


def load_top_unknown_frames(input_dir: pathlib.Path) -> polars.DataFrame:
    """Load every bcl-convert Top_Unknown_Barcodes.csv into a long barcode table.

    Expected layout: ``input_dir/<Dataset>/Reports/Top_Unknown_Barcodes.csv``;
    the dataset name is the directory that contains the ``Reports`` folder
    (``file.parents[1]``). The ``index2`` column may be missing or empty for
    single-index runs. Files below any ``Demultiplexing/`` directory are
    excluded (see ``load_stats_frames``).

    Raises:
        ValueError: If no Top_Unknown_Barcodes.csv files are found.
    """
    report_files = [
        f
        for f in input_dir.glob(f"**/Reports/{TOP_UNKNOWN_FILE_NAME}")
        if EXCLUDED_DIR not in f.parts
    ]
    if not report_files:
        raise ValueError(f"No {TOP_UNKNOWN_FILE_NAME} files found under {input_dir}.")

    frames = []
    for f in report_files:
        # Read every column as text so indexes are never parsed as numbers.
        frame = polars.read_csv(f, infer_schema=False)
        if "index2" not in frame.columns:
            frame = frame.with_columns(polars.lit("").alias("index2"))
        frames.append(
            frame.select(
                polars.col("Lane").cast(polars.Int64),
                polars.col("index").fill_null("").alias("index1"),
                polars.col("index2").fill_null(""),
                polars.col("# Reads")
                .str.replace_all(",", "")
                .cast(polars.Int64)
                .alias("Count"),
            )
            .with_columns(
                polars.when(polars.col("index2") == "")
                .then(polars.col("index1"))
                .otherwise(polars.concat_str("index1", "index2", separator=INDEX_SEP))
                .alias("Barcode"),
                polars.lit(f.parents[1].name).alias("Dataset"),
            )
            .with_columns(polars.col("Barcode").str.len_chars().alias("Barcode_len"))
            .select(
                [
                    "Lane",
                    "Barcode",
                    "Count",
                    "Dataset",
                    "Barcode_len",
                    "index1",
                    "index2",
                ]
            )
        )
    return polars.concat(frames)


def validate_table(table: polars.DataFrame) -> None:
    """Fail fast on input that violates the expected data contract.

    Every barcode must be ``index1`` or ``index1+index2`` with a non-empty
    index1, and each (Dataset, Lane) must use a single index1 length and a
    single index2 length.

    Raises:
        ValueError: If the table violates the data contract.
    """
    # index2 may be empty (a single-index barcode); index1 must not be.
    bad = table.filter(polars.col("index1").is_null() | (polars.col("index1") == ""))
    if bad.height:
        examples = bad.select("Dataset", "Barcode").unique().head(5)
        raise ValueError(
            "Found barcodes with an empty index1; every barcode must be "
            f"'index1' or 'index1{INDEX_SEP}index2'. Examples:\n{examples}"
        )
    mixed = (
        table.group_by(["Dataset", "Lane"])
        .agg(
            polars.col("index1").str.len_chars().n_unique().alias("i1_lengths"),
            polars.col("index2").str.len_chars().n_unique().alias("i2_lengths"),
        )
        .filter((polars.col("i1_lengths") > 1) | (polars.col("i2_lengths") > 1))
    )
    if mixed.height:
        raise ValueError(
            f"Expected one i1/i2 index length per (Dataset, Lane), found:\n{mixed}"
        )


def find_shared_lanes(table: polars.DataFrame) -> list[int]:
    """Return the lanes that appear in more than one dataset."""
    return (
        table.unique(subset=["Dataset", "Lane"])
        .group_by("Lane")
        .agg(polars.len().alias("Count"))
        .filter(polars.col("Count") > 1)
        .get_column("Lane")
        .to_list()
    )


def index_lengths(frame: polars.DataFrame) -> tuple[int, int]:
    """Return (i1 length, i2 length) of ``frame`` (unique by data contract)."""
    return (
        len(frame.head(1).get_column("index1").to_list()[0]),
        len(frame.head(1).get_column("index2").to_list()[0]),
    )


def has_index2(frame: polars.DataFrame) -> bool:
    """True if ``frame`` carries a non-empty second index (a dual-index dataset).

    A dataset is uniformly single-index (index2 == "") or uniformly dual-index
    (index2 non-empty) by the ``validate_table`` data contract, so one row
    decides.
    """
    return index_lengths(frame)[1] > 0


def _join_single_index(
    running: polars.DataFrame,
    other: polars.DataFrame,
    ds: str,
    count_cols: list[str],
    running_i1: int,
    other_i1: int,
) -> tuple[polars.DataFrame, bool]:
    """Join the running intersection with a step where at least one side is
    single-index, comparing on index1 only.

    Trims the longer index1 to the shorter (either side), collapses both sides
    over (index1, Lane) — summing each dataset's counts and discarding index2 —
    and inner-joins on index1 alone. Returns the new running frame (index2 set
    to "" and Barcode set to index1) and False (the result carries no index2).
    """
    count_col = f"{COUNT_COLUMN_PREFIX}{ds}"
    target_i1 = min(running_i1, other_i1)
    if running_i1 > target_i1:
        running = running.with_columns(polars.col("index1").str.slice(0, target_i1))
    if other_i1 > target_i1:
        other = other.with_columns(polars.col("index1").str.slice(0, target_i1))
    running = (
        running.group_by(["index1", "Lane"])
        .agg(
            [polars.sum(c).alias(c) for c in count_cols]
            + [polars.col("Dataset").first().alias("Dataset")]
        )
        .select(["Lane", *count_cols, "Dataset", "index1"])
    )
    other = (
        other.group_by(["index1", "Lane"])
        .agg(polars.sum(count_col).alias(count_col))
        .select(["Lane", count_col, "index1"])
    )
    return (
        running.join(other, on=["index1", "Lane"], how="inner")
        .with_columns(
            polars.col("index1").alias("Barcode"),
            polars.lit("", dtype=polars.String).alias("index2"),
        )
        .select(
            [
                "Lane",
                "Barcode",
                *count_cols,
                count_col,
                "Dataset",
                "index1",
                "index2",
            ]
        ),
        False,
    )


def intersect_lane_datasets(lane_dt: polars.DataFrame) -> polars.DataFrame:
    """Compute the barcodes undetermined in *all* of a lane's datasets.

    ``lane_dt`` must be sorted so longer barcodes come first (see ``main``);
    with ``unique(maintain_order=True)`` the first dataset is therefore the
    longest-index one (the hub), so the side that gets truncated is
    deterministic. The running intersection is inner-joined against each
    subsequent dataset, truncating and summing counts when the next dataset
    has shorter indexes.

    Datasets may be dual-index (``index1+index2``) or single-index
    (``index1`` only, i.e. index2 == ""). When every dataset in the lane is
    dual-index, the join is on (index1, index2) and the result keeps both
    indexes. As soon as a single-index dataset joins, the comparison drops to
    index1 only: the longer index1 is trimmed to the shorter (whichever side),
    index2 is discarded, and the result is reported in the single-index form.

    Each dataset's raw count is carried in its own ``Count_<dataset>`` column;
    the reported ``Count`` is their mean, rounded up to a multiple of
    ``COUNT_ROUND_TO`` with a single rounding at the end. The reported
    ``Barcode`` comes from the last-joined (shortest) dataset, i.e. it is in
    the shortest index-length space of the intersection.

    Known limitation: two datasets with equal total index length but crossed
    per-index lengths (e.g. 10+8 vs 8+10) tie in the sort and the name order
    then decides the hub, which can leave one index untruncated and make that
    lane's intersection come out empty.
    """
    lane_datasets = lane_dt["Dataset"].unique(maintain_order=True).to_list()
    by_dataset = {
        ds: lane_dt.filter(polars.col("Dataset") == ds) for ds in lane_datasets
    }

    hub = lane_datasets[0]
    running = by_dataset[hub].rename({"Count": f"{COUNT_COLUMN_PREFIX}{hub}"})
    running_has_i2 = has_index2(running)

    for ds in lane_datasets[1:]:
        other = by_dataset[ds].rename({"Count": f"{COUNT_COLUMN_PREFIX}{ds}"})
        count_cols = [c for c in running.columns if c.startswith(COUNT_COLUMN_PREFIX)]
        running_i1, running_i2 = index_lengths(running)
        other_i1, other_i2 = index_lengths(other)

        if not (running_has_i2 and has_index2(other)):
            # A single-index dataset is present (or the running intersection is
            # already single-index): compare on index1 only and drop index2.
            running, running_has_i2 = _join_single_index(
                running, other, ds, count_cols, running_i1, other_i1
            )
            if running.is_empty():
                break
            continue

        # Both dual-index: the original all-dual path, unchanged.
        if running_i1 > other_i1 or running_i2 > other_i2:
            if running_i1 > other_i1:
                running = running.with_columns(
                    polars.col("index1").str.slice(0, other_i1)
                )
            if running_i2 > other_i2:
                running = running.with_columns(
                    polars.col("index2").str.slice(0, other_i2)
                )
            # Truncation can collapse distinct full barcodes to the same
            # (index1, index2); sum each dataset's counts over the collapse.
            running = (
                running.group_by(["index1", "index2", "Lane"])
                .agg(
                    [polars.sum(c).alias(c) for c in count_cols]
                    + [
                        polars.col("Barcode").first().alias("Barcode"),
                        polars.col("Dataset").first().alias("Dataset"),
                    ]
                )
                .select(["Lane", "Barcode", *count_cols, "Dataset", "index1", "index2"])
            )

        running = (
            running.join(other, on=["index1", "index2", "Lane"], how="inner")
            .select(
                [
                    "Lane",
                    "Barcode_right",
                    *count_cols,
                    f"{COUNT_COLUMN_PREFIX}{ds}",
                    "Dataset",
                    "index1",
                    "index2",
                ]
            )
            .rename({"Barcode_right": "Barcode"})
        )
        if running.is_empty():
            break

    count_cols = [c for c in running.columns if c.startswith(COUNT_COLUMN_PREFIX)]
    mean_rounded_up = (
        polars.sum_horizontal(count_cols) / (COUNT_ROUND_TO * len(count_cols))
    ).ceil() * COUNT_ROUND_TO
    running = running.with_columns(mean_rounded_up.cast(polars.Int64).alias("Count"))
    return running.select(["Lane", "Barcode", "Count"])


def finalize(pair_frames: list[polars.DataFrame]) -> polars.DataFrame:
    """Concat the per-lane results and add the per-lane percentage column."""
    if not pair_frames:
        return polars.DataFrame(
            schema={
                "Lane": polars.Int64,
                "Barcode": polars.String,
                "Count": polars.Int64,
                "Total %": polars.Float64,
            }
        )
    return (
        polars.concat(pair_frames)
        .with_columns(
            (polars.col("Count") / polars.col("Count").sum() * 100)
            .over("Lane")
            .alias("Total %")
        )
        .sort(["Lane", "Count", "Barcode"], descending=[False, True, False])
    )


def shared_undetermined(
    table: polars.DataFrame, lanes_to_filter: list[int]
) -> polars.DataFrame:
    """Barcodes undetermined in every dataset of each lane in ``lanes_to_filter``."""
    table = (
        table.filter(polars.col("Lane").is_in(lanes_to_filter))
        # Longer barcodes first: combined with unique(maintain_order=True) this
        # makes the longest-index dataset the hub (the truncated side).
        .sort(["Lane", "Barcode_len", "Dataset"], descending=[False, True, False])
        .drop(["Barcode_len"])
    )

    validate_table(table)

    pair_frames = [
        intersect_lane_datasets(table.filter(polars.col("Lane") == lane))
        for lane in lanes_to_filter
    ]

    return finalize(pair_frames)


def restrict_lanes(lanes: list[int], requested: list[int] | None) -> list[int]:
    """Keep only the ``requested`` lanes (all when None), warning about the
    requested lanes that are not in ``lanes``."""
    if requested is None:
        return lanes
    for lane in sorted(set(requested) - set(lanes)):
        logging.warning(f"Requested lane {lane} is not eligible; skipping")
    return [lane for lane in lanes if lane in requested]


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if args.output_file is None:
        args.quiet = True
    setup_logging(args)
    logging.info("Running waif module...")

    input_dir = args.input_path

    if args.exact:
        sample_sheet = args.sample_sheet or input_dir / "SampleSheet.csv"
        lanes_to_filter = get_shared_lanes_from_samplesheet(sample_sheet)
        if not lanes_to_filter:
            raise ValueError(
                "No lane with samples from two or more projects found in the "
                "Sample Sheet; nothing to do in --exact mode."
            )
        lanes_to_filter = restrict_lanes(lanes_to_filter, args.lanes)
        if not lanes_to_filter:
            raise ValueError(
                f"None of the requested lanes {sorted(set(args.lanes))} carries "
                "samples from two or more projects; nothing to do in --exact mode."
            )
        table = fastqs_to_stats_frames(input_dir, lanes_to_filter, args.block_cache)
    else:
        if args.bclconvert_reports:
            table = load_top_unknown_frames(input_dir)
        else:
            table = load_stats_frames(input_dir)
        lanes_to_filter = restrict_lanes(find_shared_lanes(table), args.lanes)

    present_lanes = set(table.get_column("Lane").unique().to_list())
    lanes_to_filter = [lane for lane in lanes_to_filter if lane in present_lanes]
    if not lanes_to_filter:
        message = (
            "No undetermined reads found for the selected lanes"
            if args.exact
            else "No lane appears in more than one dataset"
        )
        if args.lanes is not None and not args.exact:
            message += " among the requested lanes"
        logging.info(f"{message}; nothing to report.")

    results = shared_undetermined(table, lanes_to_filter)

    if args.output_file is not None:
        results.write_csv(args.output_file)
        logging.info(f"Results saved to {args.output_file}.")
    else:
        sys.stdout.write(results.write_csv())
