"""Strainer main script."""

import csv
import dnaio
import polars
import regex
import pathlib
from collections import Counter
import argparse
import logging
from multiprocessing import Pool

from biomate.setup import setup_logging


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
    """Extract indexes information from the Sample Sheet."""
    df = []
    # Flag to track whether we are in the [Data] section of the Sample Sheet
    in_data_section = False
    # Flag to track whether the header of the [Data] section has been parsed
    header_parsed = False
    with open(sample_sheet_path, "r") as infile:
        for fields in csv.reader(infile):
            if in_data_section:
                if not header_parsed:
                    # Parse the header to find the column indices for Lane, Sample_Project, index, and index2
                    lane_idx = next(
                        (i for i, x in enumerate(fields) if x == "Lane"), None
                    )
                    sample_project_idx = next(
                        (i for i, x in enumerate(fields) if x == "Sample_Project"), None
                    )
                    index_idx = next(
                        (i for i, x in enumerate(fields) if x == "index"), None
                    )
                    index2_idx = next(
                        (i for i, x in enumerate(fields) if x == "index2"), None
                    )
                    header_parsed = True
                    if (
                        lane_idx is None
                        or sample_project_idx is None
                        or index_idx is None
                    ):
                        raise ValueError(
                            "Could not find required columns 'Lane', 'Sample_Project', and 'index' in the Sample Sheet."
                        )
                else:
                    # Parse the data lines to extract the relevant information
                    df.append(
                        {
                            "Lane": fields[lane_idx],
                            "Sample_Project": fields[sample_project_idx],
                            "index": fields[index_idx],
                            "index2": fields[index2_idx]
                            if index2_idx is not None
                            else None,
                        }
                    )
            # Trigger for entering the [Data] (v1) or [BCLConvert_Data] (v2) section of the Sample Sheet
            if fields and fields[0].strip() in ["[Data]", "[BCLConvert_Data]"]:
                in_data_section = True
                continue

    if not in_data_section:
        raise ValueError(
            "Could not find the [Data] or [BCLConvert_Data] section in the Sample Sheet."
        )
    if not header_parsed:
        raise ValueError(
            "Could not find the header line in the [Data] or [BCLConvert_Data] section of the Sample Sheet."
        )

    if not df:
        raise ValueError("No data extracted from the Sample Sheet.")

    return polars.from_dicts(df).with_columns(
        polars.col("Lane").str.zfill(3).str.pad_start(4, "L")
    )


def extract_indexes_from_undetermined_file(
    undetermined_file_path: pathlib.Path,
) -> list[dict[str, str | int]]:
    """Extract indexes information from an Undetermined FASTQ file."""
    # Extract lane information from the file name
    lane = undetermined_file_path.name.split("_")[2]
    if not len(lane) >= 3 or not regex.match(r"^L\d{3}$", lane):
        raise ValueError(
            f"Could not extract lane information from file name {undetermined_file_path.name}."
        )
    indexes = Counter()
    with dnaio.open(undetermined_file_path) as reader:
        failed_split = False
        for record in reader:
            split_name = record.name.split(" ")
            if len(split_name) < 2:
                failed_split = True
                break
            idx = split_name[1].split(":")[-1]
            if (
                "GGGG" not in idx
            ):  # Skip indexes with 4 or more consecutive Gs, which are likely to be sequencing errors
                indexes.update({idx})

    if failed_split:
        raise ValueError(
            f"Could not extract index information from file {undetermined_file_path.name}."
        )

    return [
        {"Lane": lane, "Index": idx, "Count": count}
        for idx, count in indexes.most_common(1000)
    ]


def parallel_extract_indexes_from_undetermined_files(
    undetermined_files: list[pathlib.Path], threads: int
) -> list[dict[str, str | int]]:
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
    """Search for unexpected indexes among the undetermined lanes."""
    df = []
    for lane in sample_sheet_df.get_column("Lane").unique():
        ss_lane_df = sample_sheet_df.filter(polars.col("Lane") == lane)
        und_lane_df = undetermined_df.filter(polars.col("Lane") != lane)
        # Iterate through each lane and index from the sample sheet
        for lane, proj, i1, i2 in ss_lane_df.iter_rows():
            # Iterate through each lane and most common indexes from the undetermined data frame
            for ln, i, c in und_lane_df.iter_rows():
                match_2 = None
                # Allow for up to 1 mismatch in each index
                match_1 = regex.search(r"(" + regex.escape(i1) + "){s<=1}", i)
                if i2 is not None:
                    match_2 = regex.search(r"(" + regex.escape(i2) + "){s<=1}", i)
                if match_1 and (match_2 or i2 is None):
                    df.append(
                        {
                            "Lane": ln,
                            "Index": i,
                            "Count": c,
                            "Sample_Project": proj,
                            "Lane_Project": lane,
                            "index1": i1,
                            "index2": i2,
                        }
                    )
    if not df:
        return polars.DataFrame(
            schema={
                "Lane": polars.String,
                "Index": polars.String,
                "Count": polars.Int64,
                "Sample_Project": polars.String,
                "Lane_Project": polars.String,
                "index1": polars.String,
                "index2": polars.String,
            }
        )
    return polars.from_dicts(df)


def clean_results(
    results: polars.DataFrame, sample_sheet_df: polars.DataFrame
) -> polars.DataFrame:
    """Clean the results by filtering out expected indexes."""
    return (
        results.join(
            sample_sheet_df,
            left_on=["Lane", "index1", "index2"],
            right_on=["Lane", "index", "index2"],
            how="left",
        )
        .filter(polars.col("Sample_Project_right").is_null())
        .drop(["Sample_Project_right"])
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
