"""Strainer main script."""

import argparse
import csv
import logging
import pathlib
from collections import Counter
from multiprocessing import Pool
from typing import cast

import dnaio
import polars
import regex

from biomate.setup import setup_logging

# Module-level constants
POLY_G_ARTIFACT = "GGGG"  # Sequence pattern indicating likely sequencing errors
RESULTS_SCHEMA = {
    "Lane": polars.String,
    "Undetermined_Index": polars.String,
    "Count": polars.Int64,
    "Fraction": polars.Float32,
    "Sample_Project": polars.String,
    "Lane_Project": polars.String,
    "index1": polars.String,
    "index2": polars.String,
}
REQUIRED_COLUMNS = ["Lane", "Sample_Project", "index"]
DATA_SECTION_MARKERS = ["[Data]", "[BCLConvert_Data]"]


def _find_column_indices(
    header: list[str], required_columns: list[str]
) -> dict[str, int | None]:
    """Find column indices for specified columns in a header row.

    Args:
        header: List of column names
        required_columns: List of column names to find

    Returns:
        Dictionary mapping column names to their indices (None if not found)
    """
    return {
        col: next((i for i, x in enumerate(header) if x == col), None)
        for col in required_columns
    }


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} is not a valid directory."
        )
    if not any(args.input_path.glob("**/Undetermined*R1*.fastq.gz")):
        raise argparse.ArgumentTypeError(
            f"No FASTQ files found in the input directory {args.input_path}."
        )
    if not (args.input_path / "SampleSheet.csv").is_file() and not args.sample_sheet:
        raise argparse.ArgumentTypeError(
            f"SampleSheet.csv not found in the input directory {args.input_path}."
        )
    if args.sample_sheet and not args.sample_sheet.is_file():
        raise argparse.ArgumentTypeError(
            f"Sample sheet {args.sample_sheet} does not exist."
        )
    if not args.output_path.is_dir():
        logging.info(
            f"Output directory {args.output_path} does not exist. Creating it..."
        )
        args.output_path.mkdir(parents=True, exist_ok=True)
    if args.threads < 1:
        raise argparse.ArgumentTypeError("--threads must be >= 1.")
    return args


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Evaluate indexes mixing across flowcell lanes.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Evaluate indexes mixing across flowcell lanes.",
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        required=True,
        help="Path to the input flowcell directory containing undetermined FASTQ files.",
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        required=False,
        default=pathlib.Path("./"),
        help="Path to the output directory where the results will be created (default: current directory).",
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        required=False,
        help="Path to the Sample Sheet file. If not provided, the script will look for SampleSheet.csv in the input directory.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of threads used to process the undetermined FASTQ files.",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def extract_indexes_from_sample_sheet(
    sample_sheet_path: pathlib.Path,
) -> polars.DataFrame:
    """Extract indexes information from the Sample Sheet.

    Expects a CSV file with [Data] or [BCLConvert_Data] section containing
    Lane, Sample_Project, index, and optional index2 columns.
    """
    data_rows = []
    in_data_section = False
    column_indices = {}

    with open(sample_sheet_path, "r") as infile:
        for fields in csv.reader(infile):
            # Check for section markers
            if fields and fields[0].strip() in DATA_SECTION_MARKERS:
                in_data_section = True
                continue

            if not in_data_section or not fields:
                continue

            # Parse header row to find column indices
            if not column_indices:
                column_indices = _find_column_indices(
                    fields, REQUIRED_COLUMNS + ["index2"]
                )
                if any(column_indices.get(col) is None for col in REQUIRED_COLUMNS):
                    missing = [
                        col
                        for col in REQUIRED_COLUMNS
                        if column_indices.get(col) is None
                    ]
                    raise ValueError(
                        f"Sample Sheet missing required columns: {', '.join(missing)}"
                    )
                continue

            # Parse data rows (column_indices is guaranteed to have all required columns at this point)
            lane_idx = cast(int, column_indices["Lane"])
            proj_idx = cast(int, column_indices["Sample_Project"])
            idx1_idx = cast(int, column_indices["index"])
            idx2_idx = column_indices.get("index2")

            data_rows.append(
                {
                    "Lane": fields[lane_idx],
                    "Sample_Project": fields[proj_idx],
                    "index1": fields[idx1_idx],
                    "index2": fields[idx2_idx] if idx2_idx is not None else None,
                }
            )

    if not in_data_section:
        raise ValueError(
            f"Could not find data section ({' or '.join(DATA_SECTION_MARKERS)}) in Sample Sheet."
        )

    if not data_rows:
        raise ValueError("No data rows found in Sample Sheet data section.")

    return polars.from_dicts(data_rows).with_columns(
        polars.col("Lane").str.zfill(3).str.pad_start(4, "L")
    )


def extract_indexes_from_undetermined_file(
    undetermined_file_path: pathlib.Path,
) -> list[dict[str, str | int | float]]:
    """Extract indexes information from an Undetermined FASTQ file.

    Reads the file header to extract index information, filtering out
    poly-G artifacts which indicate likely sequencing errors.
    """
    # Extract lane information from the file name
    lane = undetermined_file_path.name.split("_")[2]
    if not regex.match(r"^L\d{3}$", lane):
        raise ValueError(
            f"Invalid lane format in file name: {undetermined_file_path.name}. Expected format: *_*_LXXX_*"
        )

    total_records = 0
    indexes = Counter()
    with dnaio.open(undetermined_file_path) as reader:
        for record in reader:
            total_records += 1
            # Expected format: name parts separated by space, index in second part after last colon
            name_parts = record.name.split(" ")
            if len(name_parts) < 2:
                raise ValueError(
                    f"Invalid read name format in {undetermined_file_path.name}: {record.name}"
                )

            # Extract index from second part (format: control:barcode)
            idx = name_parts[1].split(":")[-1]

            # Skip indexes with poly-G artifacts (likely sequencing errors)
            if POLY_G_ARTIFACT not in idx:
                indexes.update({idx})

    return [
        {
            "Lane": lane,
            "Undetermined_Index": idx,
            "Count": count,
            "Fraction": count / total_records,
        }
        for idx, count in indexes.most_common(1000)
    ]


def parallel_extract_indexes_from_undetermined_files(
    undetermined_files: list[pathlib.Path], threads: int
) -> list[dict[str, str | int | float]]:
    """Extract indexes information from multiple Undetermined FASTQ files in parallel."""
    if threads > 1:
        with Pool(threads) as pool:
            results = pool.map(
                extract_indexes_from_undetermined_file, undetermined_files
            )
    else:
        results = [
            extract_indexes_from_undetermined_file(f) for f in undetermined_files
        ]

    # Flatten the list of lists of dictionaries into a single list of dictionaries
    flat_results = [item for sublist in results for item in sublist]
    return flat_results


def search_for_unexpected_indexes(
    sample_sheet_df: polars.DataFrame, undetermined_df: polars.DataFrame
) -> polars.DataFrame:
    """Search for unexpected indexes among the undetermined lanes.

    Compares undetermined indexes with expected indexes from the sample sheet,
    allowing up to 1 mismatch per index.
    """
    results = []

    for lane in sample_sheet_df.get_column("Lane").unique():
        ss_lane_df = sample_sheet_df.filter(polars.col("Lane") == lane)
        und_lane_df = undetermined_df.filter(polars.col("Lane") != lane)

        # Iterate through each sample in the sample sheet for this lane
        for ss_lane, ss_project, ss_index1, ss_index2 in ss_lane_df.iter_rows():
            # Iterate through undetermined indexes from other lanes
            for und_lane, und_index, und_count, und_fraction in und_lane_df.iter_rows():
                # Allow up to 1 mismatch in each index
                match_1 = regex.search(
                    r"(" + regex.escape(ss_index1) + "){s<=1}", und_index
                )
                match_2 = None
                if ss_index2 is not None:
                    match_2 = regex.search(
                        r"(" + regex.escape(ss_index2) + "){s<=1}", und_index
                    )

                # If both indices match (or index2 is not required), record as potential mixing
                if match_1 and (match_2 or ss_index2 is None):
                    results.append(
                        {
                            "Lane": und_lane,
                            "Undetermined_Index": und_index,
                            "Count": und_count,
                            "Fraction": und_fraction,
                            "Sample_Project": ss_project,
                            "Lane_Project": ss_lane,
                            "index1": ss_index1,
                            "index2": ss_index2,
                        }
                    )

    if not results:
        return polars.DataFrame(schema=RESULTS_SCHEMA)

    return polars.from_dicts(results)


def clean_results(
    results: polars.DataFrame, sample_sheet_df: polars.DataFrame
) -> polars.DataFrame:
    """Clean the results by filtering out expected indexes.

    Removes indexes that match the expected indexes from the sample sheet,
    keeping only truly unexpected indexes that indicate potential sample mixing.
    """
    # Join with sample sheet to identify expected indexes
    joined = results.join(
        sample_sheet_df,
        left_on=["Lane", "index1", "index2"],
        right_on=["Lane", "index1", "index2"],
        how="left",
        suffix="_expected",
    )

    # Keep only rows where the expected columns are null (i.e., unexpected indexes)
    return (
        joined.filter(polars.col("Sample_Project_expected").is_null())
        .drop("Sample_Project_expected")
        .sort(["Lane", "Count"], descending=[False, True])
    )


def aggregate_results_by_lane(results: polars.DataFrame) -> polars.DataFrame:
    """Aggregate results by lane, summing counts and fractions for each lane."""
    return (
        results.group_by(["Lane", "Undetermined_Index"])
        .agg(
            polars.col("Count").first(),
            polars.col("Fraction").first(),
            polars.col("Sample_Project"),
            polars.col("Lane_Project"),
        )
        .group_by("Lane")
        .agg(
            polars.col("Count").sum().alias("Lane_Count"),
            polars.col("Fraction").sum().alias("Lane_Fraction"),
        )
        .sort("Lane")
    )


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running strainer module...")

    # Extract the flowcell name from the input path for logging and output file naming
    basename = args.input_path.name
    logging.debug(f"Processing flowcell: '{basename}'")
    output_file = args.output_path / f"{basename}_unexpected_indexes.csv"
    logging.debug(f"Output will be saved to: {output_file}")
    if args.threads > 1:
        logging.debug(f"Using {args.threads} threads for processing.")

    # Determine the sample sheet path, either from the command line argument or by looking for SampleSheet.csv in the input directory
    sample_sheet = args.sample_sheet or args.input_path / "SampleSheet.csv"
    logging.debug(f"Using sample sheet: {sample_sheet}")

    # Extract indexes from the sample sheet
    sample_sheet_indexes = extract_indexes_from_sample_sheet(sample_sheet)

    # Recursively find all undetermined FASTQ files in the input directory
    undetermined_files = list(args.input_path.glob("**/Undetermined*R1*.fastq.gz"))
    logging.debug(f"Found {len(undetermined_files)} undetermined FASTQ files.")

    logging.info("Extracting indexes from undetermined FASTQ files...")
    undetermined_indexes = parallel_extract_indexes_from_undetermined_files(
        undetermined_files, args.threads
    )

    if not undetermined_indexes:
        logging.info("No indexes extracted from undetermined FASTQ files. Exiting.")
        return

    logging.debug("Extracting indexes from undetermined FASTQ files completed.")

    undetermined_indexes = polars.from_dicts(undetermined_indexes)

    logging.info("Searching for unexpected indexes...")
    results = search_for_unexpected_indexes(sample_sheet_indexes, undetermined_indexes)

    logging.debug("Searching for unexpected indexes completed.")

    if results.is_empty():
        logging.info("No unexpected indexes found. Exiting.")
    else:
        # Sort results by lane and count, and filter out indexes that are expected in the lane
        results = clean_results(
            results.sort(["Lane", "Count"], descending=[False, True]),
            sample_sheet_indexes,
        )
        results.write_csv(output_file)
        logging.info(f"Results saved to {output_file}.")

        aggregate_results_by_lane(results).write_csv(
            args.output_path / f"{basename}_lane_summary.csv"
        )
        logging.info(
            f"Lane summary saved to {args.output_path / f'{basename}_lane_summary.csv'}."
        )
