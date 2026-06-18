"""Strainer main script."""

import dnaio
import polars
import regex
import pathlib
from collections import Counter
import argparse
import logging


def validate_args(args) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} is not a valid directory."
        )
    if not any(args.input_path.glob("Undetermined*R1*.fastq.gz")):
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
        default=0,
        help="Number of threads used by dnaio to read/write files. Default is 0, which corresponds to a single thread.",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def extract_indexes_from_sample_sheet(
    sample_sheet_path: pathlib.Path,
) -> polars.DataFrame:
    """Extract indexes information from the Sample Sheet."""
    df = []
    in_data_section = False
    header_parsed = False
    with open(sample_sheet_path, "r") as infile:
        for line in infile:
            if in_data_section:
                if not header_parsed:
                    header = line.strip().split(",")
                    lane_idx = [i for i, x in enumerate(header) if x == "Lane"][0]
                    sample_project_idx = [
                        i for i, x in enumerate(header) if x == "Sample_Project"
                    ][0]
                    index_idx = [i for i, x in enumerate(header) if x == "index"][0]
                    index2_idx = next(
                        (i for i, x in enumerate(header) if x == "index2"), None
                    )
                    header_parsed = True
                else:
                    fields = line.strip().split(",")
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
            if line.startswith("[Data]") or line.startswith("[BCLConvert_Data]"):
                in_data_section = True
                continue

    return polars.from_dicts(df).with_columns(
        polars.col("Lane").str.zfill(3).str.pad_start(4, "L")
    )


def extract_indexes_from_undetermined_file(
    undetermined_file_path: pathlib.Path, threads: int
) -> list[tuple[str, int]]:
    """Extract indexes information from an Undetermined FASTQ file."""
    indexes = Counter()
    with dnaio.open(undetermined_file_path, open_threads=threads) as reader:
        for record in reader:
            idx = record.name.split(" ")[1].split(":")[-1]
            if (
                "GGGG" not in idx
            ):  # Skip indexes with more than 4 consecutive Gs, which are likely to be sequencing errors
                indexes.update({idx})

    return indexes.most_common(1000)


def search_for_unexpected_indexes(
    sample_sheet_df: polars.DataFrame, undetermined_df: polars.DataFrame
) -> polars.DataFrame:
    """Search for unexpected indexes among the undetermined lanes."""
    df = list()
    for lane, proj, i1, i2 in sample_sheet_df.iter_rows():
        for l, i, c in undetermined_df.iter_rows():
            match_1 = regex.search(r"(" + i1 + "){s<=1}", i)
            if i2 is not None:
                match_2 = regex.search(r"(" + i2 + "){s<=1}", i)
            if match_1 and (match_2 or i2 is None):
                if lane != l:
                    df.append(
                        {
                            "lane": l,
                            "Index": i,
                            "Count": c,
                            "Sample_Project": proj,
                            "Lane_Project": lane,
                            "index1": i1,
                            "index2": i2,
                        }
                    )

    return polars.from_dicts(df)


def main(args: argparse.Namespace) -> None:
    """Main function."""
    # setup_logging(args)
    logging.info("Running strainer module...")

    basename = args.input_path.name

    sample_sheet = args.sample_sheet or args.input_path.joinpath("SampleSheet.csv")

    sample_sheet_indexes = extract_indexes_from_sample_sheet(sample_sheet)

    undetermined_files = list(args.input_path.glob("Undetermined*R1*.fastq.gz"))
    undetermined_indexes = []
    for un_file in undetermined_files:
        lane = un_file.name.split("_")[2]
        undetermined_indexes += [
            {"Lane": lane, "Index": idx, "Count": count}
            for idx, count in extract_indexes_from_undetermined_file(
                un_file, args.threads
            )
        ]
    undetermined_indexes = polars.from_dicts(undetermined_indexes)

    search_for_unexpected_indexes(sample_sheet_indexes, undetermined_indexes).write_csv(
        args.output_path.joinpath(f"{basename}_unexpected_indexes.csv")
    )
