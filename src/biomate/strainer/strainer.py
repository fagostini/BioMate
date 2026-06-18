"""Strainer main script."""

import dnaio
import polars
import regex
import pathlib
from collections import Counter
import argparse
import logging
from multiprocessing import Pool

from biomate.setup import setup_logging


def validate_args(args) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} is not a valid directory."
        )
    if not any(args.input_path.glob("**/Undetermined*R1*.fastq.gz")):
        raise argparse.ArgumentTypeError(
            f"No FASTQ files found in the input directory {args.input_path}."
        )
    if (
        not args.input_path.joinpath("SampleSheet.csv").is_file()
        and not args.sample_sheet
    ):
        raise argparse.ArgumentTypeError(
            f"SampleSheet.csv not found in the input directory {args.input_path}."
        )
    if args.sample_sheet and not args.sample_sheet.is_file():
        raise argparse.ArgumentTypeError(
            f"Sample sheet {args.sample_sheet} does not exist."
        )
    if not args.output_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Output path {args.output_path} does not exist or is not a valid directory."
        )
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
        for line in infile:
            if in_data_section:
                if not header_parsed:
                    # Parse the header to find the column indices for Lane, Sample_Project, index, and index2
                    header = line.strip().split(",")
                    lane_idx = [i for i, x in enumerate(header) if x == "Lane"][0]
                    sample_project_idx = [
                        i for i, x in enumerate(header) if x == "Sample_Project"
                    ][0]
                    index_idx = [i for i, x in enumerate(header) if x == "index"][0]
                    index2_idx = [i for i, x in enumerate(header) if x == "index2"][0]
                    header_parsed = True
                else:
                    # Parse the data lines to extract the relevant information
                    fields = line.strip().split(",")
                    df.append(
                        {
                            "Lane": fields[lane_idx],
                            "Sample_Project": fields[sample_project_idx],
                            "index": fields[index_idx],
                            "index2": fields[index2_idx],
                        }
                        if index2_idx
                        else {
                            "Lane": fields[lane_idx],
                            "Sample_Project": fields[sample_project_idx],
                            "index": fields[index_idx],
                        }
                    )
            # Trigger for entering the [Data] (v1) or [BCLConvert_Data] (v2) section of the Sample Sheet
            if line.startswith("[Data]") or line.startswith("[BCLConvert_Data]"):
                in_data_section = True
                continue

    return polars.from_dicts(df).with_columns(
        polars.col("Lane").str.zfill(3).str.pad_start(4, "L")
    )


def extract_indexes_from_undetermined_file(
    undetermined_file_path: pathlib.Path,
) -> list[dict[str, str | int]]:
    """Extract indexes information from an Undetermined FASTQ file."""
    # Extract lane information from the file name
    lane = undetermined_file_path.name.split("_")[2]
    indexes = Counter()
    with dnaio.open(undetermined_file_path) as reader:
        for record in reader:
            idx = record.name.split(" ")[1].split(":")[-1]
            # Filter out indexes containing "GGGG" as they are likely to be sequencing errors
            if "GGGG" not in indexes:
                indexes.update({idx})

    return [
        {"Lane": lane, "Index": idx, "Count": count}
        for idx, count in indexes.most_common(1000)
    ]


def search_for_unexpected_indexes(
    sample_sheet_df: polars.DataFrame, undetermined_df: polars.DataFrame
) -> polars.DataFrame:
    """Search for unexpected indexes among the undetermined lanes."""
    df = list()
    # Iterate through each lane and index from the sample sheet
    for lane, proj, i1, i2 in sample_sheet_df.iter_rows():
        # Iterate through each lane and most common indexes from the undetermined data frame
        for ln, i, c in undetermined_df.iter_rows():
            # Allow for up to 1 mismatch in both index 1 and index 2
            if regex.search(r"(" + i1 + "){s<=1}", i) and regex.search(
                r"(" + i2 + "){s<=1}", i
            ):
                if lane != ln:
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
                "Lane": polars.Utf8,
                "Index": polars.Utf8,
                "Count": polars.Int64,
                "Sample_Project": polars.Utf8,
                "Lane_Project": polars.Utf8,
                "index1": polars.Utf8,
                "index2": polars.Utf8,
            }
        )
    return polars.from_dicts(df)


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running strainer module...")

    # Extract the flowcell name from the input path for logging and output file naming
    basename = args.input_path.name
    logging.debug(f"Processing flowcell: '{basename}'")
    output_file = args.output_path.joinpath(f"{basename}_unexpected_indexes.csv")
    logging.debug(f"Output will be saved to: {output_file}")
    if args.threads > 1:
        logging.debug(f"Using {args.threads} threads for processing.")

    # Determine the sample sheet path, either from the command line argument or by looking for SampleSheet.csv in the input directory
    sample_sheet = args.sample_sheet or args.input_path.joinpath("SampleSheet.csv")
    logging.debug(f"Using sample sheet: {sample_sheet}")

    # Extract indexes from the sample sheet
    sample_sheet_indexes = extract_indexes_from_sample_sheet(sample_sheet)

    # Recursively find all undetermined FASTQ files in the input directory
    undetermined_files = list(args.input_path.glob("**/Undetermined*R1*.fastq.gz"))
    logging.debug(f"Found {len(undetermined_files)} undetermined FASTQ files.")

    logging.info("Extracting indexes from undetermined FASTQ files...")
    undetermined_indexes = []
    # Use multiprocessing to extract indexes from undetermined FASTQ files in parallel
    with Pool(args.threads) as pool:
        undetermined_indexes = pool.map(
            extract_indexes_from_undetermined_file, undetermined_files
        )
    # Flatten the list of lists of dictionaries into a single list of dictionaries and convert it to a Polars DataFrame
    undetermined_indexes = polars.from_dicts(
        [item for sublist in undetermined_indexes for item in sublist]
    )

    logging.info("Searching for unexpected indexes...")
    results = search_for_unexpected_indexes(sample_sheet_indexes, undetermined_indexes)

    if results.is_empty():
        logging.info("No unexpected indexes found. Exiting.")
        results.write_csv(output_file)
    else:
        # Sort results by lane and count, and filter out indexes that are expected in the lane
        results.sort(["Lane", "Count"], descending=[False, True]).join(
            sample_sheet_indexes,
            left_on=["Lane", "index1", "index2"],
            right_on=["Lane", "index", "index2"],
            how="left",
        ).filter(polars.col("Sample_Project_right").is_null()).drop(
            ["Sample_Project_right"]
        ).write_csv(output_file)
