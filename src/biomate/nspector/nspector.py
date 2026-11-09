"""Nspector main script."""

import argparse
import logging
import pathlib

import altair as alt
import polars

from biomate.nspector.report import write_multiqc, write_report
from biomate.nspector.scan import (
    FileScan,
    extract_tile_coordinates,  # noqa: F401  (re-exported for callers/tests)
    scan_files,
)
from biomate.nspector.tables import build_tables, reads_frame, tile_layout
from biomate.setup import setup_logging


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    for input_file in args.input:
        if not input_file.exists():
            raise argparse.ArgumentTypeError(f"Input file {input_file} does not exist.")
        elif not input_file.is_file():
            raise argparse.ArgumentTypeError(f"Input path {input_file} is not a file.")

    if getattr(args, "spatial_bins", 1) < 1:
        raise argparse.ArgumentTypeError("--spatial-bins must be at least 1.")

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
        description="Inspect sequencing data for N content and their distribution "
        "across lanes, samples, flowcell tiles and cycles. Writes an interactive "
        "HTML report, tidy Parquet tables and MultiQC custom content.",
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
    parser.add_argument(
        "--engine",
        type=str,
        default="python",
        choices=["python", "rust"],
        help="Scan engine: pure Python (default) or the native biomate_kernels "
        "Rust scan (falls back to Python if biomate_kernels is not installed).",
    )
    parser.add_argument(
        "--spatial-bins",
        type=int,
        default=50,
        help="Grid size (bins per axis) for the per-tile spatial N maps (default: 50).",
    )
    parser.add_argument(
        "--png",
        action="store_true",
        help="Also write static PNG charts per input file (slow for many tiles).",
    )
    parser.add_argument(
        "--no-html",
        action="store_true",
        help="Do not write the interactive HTML report.",
    )
    parser.add_argument(
        "--no-multiqc",
        action="store_true",
        help="Do not write MultiQC custom-content files.",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


# ---------------------------------------------------------------------------
# Static PNG charts (opt-in with --png)
# ---------------------------------------------------------------------------


def exploded_frame(scan: FileScan) -> polars.DataFrame:
    """One row per N base: tile, x, y and 1-based cycle."""
    return (
        reads_frame(scan)
        .select("tile", "x", "y", "cycles")
        .explode("cycles")
        .with_columns(polars.col("cycles").cast(polars.UInt16))
    )


def plot_n_per_sequence(data: polars.DataFrame, output_path: pathlib.Path) -> None:
    """Plot the distribution of N counts per sequence."""
    chart = (
        data.group_by(["tile", "x", "y"])
        .agg(polars.len().alias("n_count"))
        .group_by("n_count")
        .agg(polars.len().alias("n_sequences"))
        .sort("n_count")
        .plot.bar(
            x=alt.X("n_count:O", title="Number of N bases"),
            y=alt.Y("n_sequences", title="Number of Sequences").scale(type="symlog"),
        )
        .properties(
            width=800, height=400, title="Distribution of N counts per sequence"
        )
    )
    chart.save(output_path)


def plot_n_per_cycle_and_tile(
    data: polars.DataFrame, output_path: pathlib.Path, max_cycle: int | None = None
) -> None:
    """Plot the distribution of N-containing sequences across cycles and tiles."""
    max_cycle = max_cycle or int(data.get_column("cycles").max())
    grid = polars.DataFrame({"tile": data.get_column("tile").unique().sort()}).join(
        polars.DataFrame(
            {"cycles": polars.arange(1, max_cycle + 1, eager=True).cast(polars.UInt16)}
        ),
        how="cross",
    )
    source = (
        grid.join(
            data.group_by(["cycles", "tile"]).agg(polars.len().alias("n_sequences")),
            on=["cycles", "tile"],
            how="left",
        )
        .with_columns(polars.col("n_sequences").fill_null(0))
        .with_columns(*tile_layout(polars.col("tile")))
        # fixed row order so the surface/swath/segment category order (and
        # hence the rendered facet layout) is deterministic across runs
        .sort(["surface", "swath", "segment", "cycles"])
    )

    tick_step = 5 if max_cycle <= 160 else 10
    chart = (
        source.plot.rect(
            x=alt.X(
                "cycles:O",
                title="Cycle",
                axis=alt.Axis(values=list(range(tick_step, max_cycle + 1, tick_step))),
            ),
            y=alt.Y("segment:O", title="Tile"),
            color=alt.Color("n_sequences:Q", title="Number of Sequences")
            .scale(scheme="viridis")
            .bin(maxbins=20),
            stroke=alt.value("white"),
            column=alt.Column("surface:N").title("Surface"),
            row=alt.Row("swath:N").title("Swath").header(labelAngle=0),
        )
        .configure_axis(labelFontSize=12)
        .properties(
            width=max_cycle * 6,
            height=source.get_column("segment").n_unique() * 12,
            title="Distribution of N-containing sequences across cycles and tiles",
        )
    )
    chart.save(output_path)


def _tile_chart(source: polars.DataFrame) -> alt.VConcatChart:
    """Build the per-tile composite chart (top histogram, points, right histogram).

    ``source`` must already be sorted by (x, y) so that the order in which
    overlapping semi-transparent bars are drawn is deterministic.
    """
    base = alt.Chart(source)
    base_bar = base.mark_bar(opacity=0.3, binSpacing=0)

    xscale = alt.Scale(
        domainMin=int(source.get_column("x").min()),
        domainMax=int(source.get_column("x").max()),
    )
    yscale = alt.Scale(
        domainMin=int(source.get_column("y").min()),
        domainMax=int(source.get_column("y").max()),
    )

    points = (
        base.mark_circle()
        .encode(alt.X("x").scale(xscale), alt.Y("y").scale(yscale))
        .properties(width=500, height=500)
    )
    top_hist = base_bar.encode(
        # when using bins, the axis scale is set through the bin extent
        alt.X("x:Q").bin(maxbins=200, extent=xscale.domain).stack(None).title(""),
        alt.Y("count()").stack(None).title(""),
    ).properties(width=500, height=60)
    right_hist = base_bar.encode(
        alt.Y("y:Q").bin(maxbins=200, extent=yscale.domain).stack(None).title(""),
        alt.X("count()").stack(None).title(""),
    ).properties(width=60, height=500)

    return top_hist & (points | right_hist)


def plot_individual_tiles(
    data: polars.DataFrame, output_path: pathlib.Path, stem: str
) -> None:
    """Plot the distribution of N-containing sequences for individual tiles.

    Saves one composite chart per tile as
    ``individual_tiles/{stem}_tile_{tile}.png``. Tiles are processed in sorted
    order and each per-tile source is sorted by (x, y), which makes the
    rendered images deterministic across runs.
    """
    output_path.joinpath("individual_tiles").mkdir(parents=True, exist_ok=True)

    for tile in sorted(data.get_column("tile").unique().to_list()):
        source = (
            data.filter(polars.col("tile") == tile)
            .group_by(["x", "y"])
            .agg(polars.len().alias("z"))
            .sort(["x", "y"])
        )
        _tile_chart(source).save(
            output_path / "individual_tiles" / f"{stem}_tile_{tile}.png"
        )


def write_pngs(scan: FileScan, output: pathlib.Path) -> None:
    """Write the static PNG charts for one input file."""
    data = exploded_frame(scan)
    plot_n_per_sequence(data, output / f"{scan.stem}_n_count_distribution.png")
    plot_n_per_cycle_and_tile(
        data,
        output / f"{scan.stem}_n_distribution_per_cycle_and_tile.png",
        max_cycle=scan.max_len,
    )
    plot_individual_tiles(data, output, scan.stem)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def usable_scans(scans: list[FileScan]) -> list[FileScan]:
    """Drop files with read-name format errors or no reads, with a warning."""
    kept = []
    for scan in scans:
        if scan.format_error is not None:
            logging.warning(
                f"Skipping file {scan.path.name} due to read name format issues "
                f"(record {scan.format_error})."
            )
        elif scan.total_reads == 0:
            logging.warning(f"No reads found in {scan.path}. Skipping analysis.")
        else:
            if not scan.tiles:
                logging.warning(f"No N bases found in {scan.path}.")
            kept.append(scan)
    return kept


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running nspector module...")

    scans = usable_scans(scan_files(args.input, args.engine, args.threads))
    if not scans:
        logging.warning("No usable input files; nothing to report.")
        return

    for scan in scans:
        logging.debug(
            f"{scan.path.name}: {len(scan.tiles)} of {scan.total_reads} reads have N bases."
        )
        if scan.tiles:
            exploded_frame(scan).write_parquet(
                args.output / f"{scan.stem}_n_sequences.parquet"
            )

    tables = build_tables(scans, spatial_bins=args.spatial_bins)
    tables.write_parquet(args.output)

    if not getattr(args, "no_html", False):
        path = write_report(scans, tables, args.output)
        logging.info(f"Interactive report written to {path}")
    if not getattr(args, "no_multiqc", False):
        write_multiqc(tables, args.output)
    if getattr(args, "png", False):
        for scan in scans:
            if scan.tiles:
                logging.debug(f"Creating PNG charts for file: {scan.path}")
                write_pngs(scan, args.output)

    logging.info("Nspector module finished.")
