"""Nspector main script."""

import argparse
import logging
import pathlib

import dnaio
import polars
import altair as alt

from biomate.setup import setup_logging


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate the command line arguments.

    """
    for input_file in args.input:
        if not input_file.exists():
            raise argparse.ArgumentTypeError(f"Input file {input_file} does not exist.")
        elif not input_file.is_file():
            raise argparse.ArgumentTypeError(f"Input path {input_file} is not a file.")

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
        tiles = []
        x_coords = []
        y_coords = []
        cycles = []
        skip_file = False
        with dnaio.open(input_file) as reader:
            for record in reader:
                matches = [i for i, c in enumerate(record.sequence) if c == "N"]
                if matches:
                    # Illumina read name format: instrument:run:flowcell:lane:tile:x:y
                    parts = record.name.split(" ")[0].split(":")
                    if len(parts) < 7:
                        skip_file = True
                        break
                    tiles.append(int(parts[4]))
                    x_coords.append(int(parts[5]))
                    y_coords.append(int(parts[6]))
                    cycles.append(matches)

        if skip_file:
            logging.warning(
                f"Skipping file {input_file.name} due to read name format issues."
            )
            continue

        if not tiles:
            logging.warning(f"No N bases found in {input_file}. Skipping analysis.")
            continue

        logging.debug(f"Found {len(tiles)} reads with N bases in {input_file}.")

        data = (
            polars.DataFrame(
                {
                    "tile": tiles,
                    "x": x_coords,
                    "y": y_coords,
                    "cycles": cycles,
                }
            )
            .explode("cycles")
            .with_columns(
                polars.col("tile").cast(polars.UInt16),
                polars.col("x").cast(polars.UInt32),
                polars.col("y").cast(polars.UInt32),
                polars.col("cycles").cast(polars.UInt16),
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
                .alias("n_sequences"),
                polars.col("cycles_right").alias("cycles").cast(polars.UInt16),
                polars.col("tile_right").alias("tile").cast(polars.UInt16),
            )
            .select(["cycles", "tile", "n_sequences"])
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
