"""Nspector main script."""

import argparse
import logging
import pathlib

import dnaio
import polars
import altair as alt

from biomate.setup import setup_logging


def extract_tile_coordinates(read_name: str) -> tuple[int, int, int] | None:
    """Extract tile and x, y coordinates from Illumina read name.

    Illumina read name format: instrument:run:flowcell:lane:tile:x:y
    Returns (tile, x, y) or None if name format is invalid.
    """
    try:
        # Get the first part before any spaces (some formats have additional info)
        base_name = read_name.split(" ")[0]
        parts = base_name.split(":")
        if len(parts) < 7:
            return None
        return (int(parts[4]), int(parts[5]), int(parts[6]))
    except (IndexError, ValueError):
        return None


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    for input_file in args.input:
        if not input_file.exists():
            raise argparse.ArgumentTypeError(f"Input file {input_file} does not exist.")
        elif not input_file.is_file():
            raise argparse.ArgumentTypeError(f"Input path {input_file} is not a file.")

    if args.output.exists() and not args.output.is_dir():
        raise argparse.ArgumentTypeError(
            f"Output path {args.output} must be a directory."
        )
    elif not args.output.exists():
        logging.info(f"Output directory {args.output} does not exist. Creating it...")
        args.output.mkdir(parents=True, exist_ok=True)
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


def generate_cycles_filler(data: polars.DataFrame, cycles: int) -> polars.DataFrame:
    """Generate a DataFrame to fill in missing cycle counts."""
    return (
        polars.DataFrame(
            {
                "tile": data.group_by(["cycles", "tile"])
                .agg(polars.len().alias("n_sequences"))
                .get_column("tile")
                .unique()
                .to_list()
            }
        )
        .with_columns(
            [
                polars.lit(list(range(1, cycles + 1))).alias("cycles"),
                polars.lit(0).alias("n_filler"),
            ]
        )
        .explode("cycles")
    )


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running nspector module...")

    for input_file in args.input:
        logging.debug(f"Processing file: {input_file}")
        records_with_n = []
        has_format_error = False

        with dnaio.open(input_file) as reader:
            for record in reader:
                n_positions = [i for i, c in enumerate(record.sequence) if c == "N"]
                if n_positions:
                    coords = extract_tile_coordinates(record.name)
                    if coords is None:
                        logging.warning(
                            f"Skipping file {input_file.name} due to read name format issues."
                        )
                        has_format_error = True
                        break
                    tile, x, y = coords
                    records_with_n.append(
                        {
                            "tile": tile,
                            "x": x,
                            "y": y,
                            "cycles": n_positions,
                        }
                    )

        # Skip if we encountered format errors
        if has_format_error:
            continue

        if not records_with_n:
            logging.warning(f"No N bases found in {input_file}. Skipping analysis.")
            continue

        logging.debug(
            f"Found {len(records_with_n)} reads with N bases in {input_file}."
        )

        data = (
            polars.DataFrame(records_with_n)
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

        logging.debug(f"Creating charts for file: {input_file}")

        chart = (
            data.group_by(["tile", "x", "y"])
            .agg(polars.len().alias("n_count"))
            .group_by("n_count")
            .agg(polars.len().alias("n_sequences"))
            .plot.bar(
                x=alt.X("n_count", title="Number of N bases"),
                y=alt.Y("n_sequences", title="Number of Sequences").scale(
                    type="symlog"
                ),
            )
            .properties(
                width=1024, height=480, title="Distribution of N counts per sequence"
            )
        )
        chart.save(args.output / f"{input_file.stem}_n_count_distribution.png")

        filler = generate_cycles_filler(
            data, int(max(data.get_column("cycles").to_list()))
        )
        chart = (
            data.group_by(["cycles", "tile"])
            .agg(polars.len().alias("n_sequences"))
            .join(
                filler,
                on=["cycles", "tile"],
                how="outer",
            )
            .with_columns(
                polars.when(polars.col("n_sequences").is_not_null())
                .then(polars.col("n_sequences"))
                .otherwise(polars.col("n_filler"))
                .alias("n_sequences_filled")
            )
            .select(
                [
                    "cycles",
                    "tile",
                    polars.col("n_sequences_filled").alias("n_sequences"),
                ]
            )
            .plot.line(
                x=alt.X("cycles", title="Cycle"),
                y=alt.Y("n_sequences", title="Number of Sequences").scale(
                    type="symlog"
                ),
                color=alt.Color("tile:N", title="Tile"),
            )
            .properties(
                width=800,
                height=400,
                title="Distribution of N-containing sequences across cycles and tiles",
            )
        )
        chart.save(args.output / f"{input_file.stem}_n_cycle_tile_distribution.png")

    logging.info("Nspector module finished.")
