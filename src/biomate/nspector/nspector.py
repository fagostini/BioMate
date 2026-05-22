"""Nspector main script."""

import argparse
import logging
import pathlib

import regex
import dnaio
import polars

from biomate.setup import setup_logging


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.

    """
    if not args.input.exist():
        raise argparse.ArgumentTypeError(f"Input file {args.input} does not exist.")
    elif not args.input.is_file():
        raise argparse.ArgumentTypeError(f"Input path {args.input} is not a file.")

    if not args.output.exists():
        logging.info(f"Output directory {args.output} does not exist. Creating it...")
        args.output.mkdir(parents=True, exist_ok=True)
    else:
        if not args.output.is_dir():
            raise argparse.ArgumentTypeError(
                f"Output path {args.output} must be a directory."
            )
        else:
            logging.warning(
                f"Output directory {args.output} already exists. Results may be overwritten."
            )
    return args


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Inspect sequencing data for N content and their distribution across the flowcell tiles.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Inspect sequencing data for N content and their distribution across the flowcell tiles.",
    )
    parser.add_argument(
        "--input",
        type=pathlib.Path,
        required=True,
        nargs="+",
        help="Path to the input FASTQ file(s) to inspect.",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=pathlib.Path("./output"),
        help="Path to the output directory where the results will be saved (default: ./output).",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if not args.output:
        args.quiet = False
    setup_logging(args)
    logging.info("Running nspector module...")

    dicts_list = []

    with dnaio.open(args.input) as reader:
        for record in reader:
            matches = list(regex.finditer("N", record.sequence))
            if matches:
                dicts_list.append(
                    dict(
                        zip(
                            ["tile", "x", "y", "cycles"],
                            record.name.split(" ")[0].split(":")[4:7]
                            + [[m.start() for m in matches]],
                        )
                    )
                )
    data = (
        polars.from_dicts(dicts_list)
        .explode("cycles")
        .with_columns(
            [
                polars.col("tile").cast(polars.UInt16),
                polars.col("x").cast(polars.UInt32),
                polars.col("y").cast(polars.UInt32),
                polars.col("cycles").cast(polars.UInt16),
            ]
        )
    )
