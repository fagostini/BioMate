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
        "--threads",
        type=int,
        default=0,
        help="Number of threads to use for reading FASTQ files (default: 0, which means no extra threads apart from the main thread).",
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


def plot_n_per_sequence(data: polars.DataFrame, output_path: pathlib.Path) -> None:
    """Plot the distribution of N counts per sequence."""
    chart = (
        data.group_by(["tile", "x", "y"])
        .agg(polars.len().alias("n_count"))
        .group_by("n_count")
        .agg(polars.len().alias("n_sequences"))
        .plot.bar(
            x=alt.X("n_count", title="Number of N bases"),
            y=alt.Y("n_sequences", title="Number of Sequences").scale(type="symlog"),
        )
        .properties(
            width=1600, height=400, title="Distribution of N counts per sequence"
        )
    )
    chart.save(output_path)


def plot_n_per_cycle_and_tile(
    data: polars.DataFrame, output_path: pathlib.Path
) -> None:
    """Plot the distribution of N-containing sequences across cycles and tiles."""
    filler = generate_cycles_filler(data, int(max(data.get_column("cycles").to_list())))

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
                "cycles_right",
                "tile_right",
                polars.col("n_sequences_filled").alias("n_sequences"),
            ]
        )
        .rename({"cycles_right": "cycles", "tile_right": "tile"})
        .drop_nulls()
        .with_columns(
            polars.col("tile")
            .cast(polars.Utf8)
            .str.slice(0, 1)
            .cast(polars.Int16)
            .alias("surface"),
            polars.col("tile")
            .cast(polars.Utf8)
            .str.slice(1, 1)
            .cast(polars.Int16)
            .alias("swath"),
            polars.col("tile")
            .cast(polars.Utf8)
            .str.slice(2, 2)
            .cast(polars.Int16)
            .alias("segment"),
        )
    )

    height = chart.get_column("segment").n_unique() * 12
    width = chart.get_column("cycles").n_unique() * 6

    chart = (
        chart.plot.rect(
            x=alt.X(
                "cycles:O", title="Cycle", axis=alt.Axis(values=list(range(0, 301, 5)))
            ),
            y=alt.Y("segment:O", title="Tile"),
            color=alt.Color("n_sequences:Q", title="Number of Sequences")
            .scale(scheme="viridis")
            .bin(maxbins=20),
            stroke=alt.value("white"),
            column=alt.Column("surface:N").title("Surface"),
            row=alt.Row("swath:N").title("Swath").header(labelAngle=0),
        )
        .configure_axis(
            labelFontSize=12,
        )
        .properties(
            width=width,
            height=height,
            title="Distribution of N-containing sequences across cycles and tiles",
        )
    )
    chart.save(output_path)


def plot_individual_tiles(data: polars.DataFrame, output_path: pathlib.Path) -> None:
    """Plot the distribution of N-containing sequences for individual tiles."""
    output_path.joinpath("individual_tiles").mkdir(parents=True, exist_ok=True)

    for tile in data.get_column("tile").unique().to_list():
        source = (
            data.filter(polars.col("tile") == tile)
            .group_by(["x", "y"])
            .agg(polars.len().alias("z"))
        )

        base = alt.Chart(source)
        base_bar = base.mark_bar(opacity=0.3, binSpacing=0)

        xscale = alt.Scale(
            domainMin=min(source.get_column("x").to_list()),
            domainMax=max(source.get_column("x").to_list()),
        )
        yscale = alt.Scale(
            domainMin=min(source.get_column("y").to_list()),
            domainMax=max(source.get_column("y").to_list()),
        )

        points = (
            base.mark_circle()
            .encode(
                alt.X("x").scale(xscale),
                alt.Y("y").scale(yscale),
                # color="Species",
            )
            .properties(width=500, height=500)
        )

        top_hist = base_bar.encode(
            alt.X("x:Q")
            # when using bins, the axis scale is set through
            # the bin extent, so we do not specify the scale here
            # (which would be ignored anyway)
            .bin(maxbins=200, extent=xscale.domain)
            .stack(None)
            .title(""),
            alt.Y("count()").stack(None).title(""),
            # alt.Color("Species:N"),
        ).properties(width=500, height=60)

        right_hist = base_bar.encode(
            alt.Y("y:Q").bin(maxbins=200, extent=yscale.domain).stack(None).title(""),
            alt.X("count()").stack(None).title(""),
            # alt.Color("Species:N"),
        ).properties(width=60, height=500)

        chart = top_hist & (points | right_hist)

        chart.save(output_path / "individual_tiles" / f"tile_{tile}_n_distribution.png")


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running nspector module...")

    for input_file in args.input:
        logging.debug(f"Processing file: {input_file}")
        records_with_n = []
        has_format_error = False

        with dnaio.open(input_file, open_threads=args.threads) as reader:
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
        data.write_parquet(args.output / f"{input_file.stem}_n_sequences.parquet")

        logging.debug(f"Creating charts for file: {input_file}")

        plot_n_per_sequence(
            data, args.output / f"{input_file.stem}_n_count_distribution.png"
        )

        plot_n_per_cycle_and_tile(
            data,
            args.output / f"{input_file.stem}_n_distribution_per_cycle_and_tile.png",
        )

        plot_individual_tiles(data, args.output)

    logging.info("Nspector module finished.")
