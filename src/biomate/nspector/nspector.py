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
from biomate.nspector.tables import build_tables, tile_layout
from biomate.setup import setup_logging


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    for input_file in args.input:
        if not input_file.exists():
            raise argparse.ArgumentTypeError(f"Input file {input_file} does not exist.")
        elif not input_file.is_file():
            raise argparse.ArgumentTypeError(f"Input path {input_file} is not a file.")

    if getattr(args, "spatial_bin_size", 1) < 1:
        raise argparse.ArgumentTypeError("--spatial-bin-size must be at least 1.")

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
        default="auto",
        choices=["auto", "python", "rust"],
        help="Scan engine: auto (default; Rust when biomate_kernels is installed, "
        "else Python), python, or rust (falls back to Python if unavailable).",
    )
    parser.add_argument(
        "--spatial-bin-size",
        type=int,
        default=1000,
        help="Width of the spatial bins in read-name X/Y coordinate units (default: 1000).",
    )
    parser.add_argument(
        "--raw",
        action="store_true",
        help="Also write every N-containing read (lane, tile, x, y, cycles) to "
        "<stem>_n_sequences/ Parquet parts. Large; uses the Python engine.",
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


def _scan_frame(scan: FileScan, field: str, columns: list[str]) -> polars.DataFrame:
    return polars.DataFrame(getattr(scan, field), schema=columns, orient="row")


def plot_n_per_sequence(scan: FileScan, output_path: pathlib.Path) -> None:
    """Plot the distribution of N counts per sequence."""
    chart = (
        _scan_frame(scan, "n_per_read", ["lane", "n_count", "n_sequences"])
        .group_by("n_count")
        .agg(polars.col("n_sequences").sum())
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


def plot_n_per_cycle_and_tile(scan: FileScan, output_path: pathlib.Path) -> None:
    """Plot the distribution of N-containing sequences across cycles and tiles."""
    max_cycle = scan.max_len
    tiles = _scan_frame(scan, "tiles", ["lane", "tile", "total", "n", "n_bases"])
    grid = tiles.select(polars.col("tile").unique().sort()).join(
        polars.DataFrame({"cycles": polars.arange(1, max_cycle + 1, eager=True)}),
        how="cross",
    )
    source = (
        grid.join(
            _scan_frame(scan, "tile_cycle", ["lane", "tile", "cycles", "n_sequences"])
            .group_by("tile", "cycles")
            .agg(polars.col("n_sequences").sum()),
            on=["tile", "cycles"],
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


def _tile_chart(source: polars.DataFrame, bin_size: int) -> alt.VConcatChart:
    """Build the per-tile composite chart: binned density plus marginal histograms.

    ``source`` must already be sorted by (x_bin, y_bin) so the rendering is
    deterministic.
    """
    base = alt.Chart(source).transform_calculate(
        x0=f"datum.x_bin * {bin_size}",
        x1=f"(datum.x_bin + 1) * {bin_size}",
        y0=f"datum.y_bin * {bin_size}",
        y1=f"(datum.y_bin + 1) * {bin_size}",
    )
    density = (
        base.mark_rect()
        .encode(
            alt.X("x0:Q", title="X").scale(zero=False),
            alt.X2("x1:Q"),
            alt.Y("y0:Q", title="Y").scale(zero=False),
            alt.Y2("y1:Q"),
            alt.Color("n:Q", title="Reads with N").scale(scheme="viridis"),
        )
        .properties(width=500, height=500)
    )
    top_hist = (
        base.mark_bar()
        .encode(alt.X("x_bin:O", axis=None), alt.Y("sum(n):Q", title=""))
        .properties(width=500, height=60)
    )
    right_hist = (
        base.mark_bar()
        .encode(alt.Y("y_bin:O", axis=None), alt.X("sum(n):Q", title=""))
        .properties(width=60, height=500)
    )
    return top_hist & (density | right_hist)


def plot_individual_tiles(
    scan: FileScan, output_path: pathlib.Path, bin_size: int
) -> None:
    """Plot the spatial distribution of N-containing reads for each tile.

    Saves one composite chart per tile as
    ``individual_tiles/{stem}_tile_{tile}.png``, in sorted tile order.
    """
    output_path.joinpath("individual_tiles").mkdir(parents=True, exist_ok=True)
    spatial = _scan_frame(scan, "spatial", ["lane", "tile", "x_bin", "y_bin", "n"])
    for tile in sorted(spatial.get_column("tile").unique().to_list()):
        source = (
            spatial.filter(polars.col("tile") == tile)
            .group_by("x_bin", "y_bin")
            .agg(polars.col("n").sum())
            .sort("x_bin", "y_bin")
        )
        _tile_chart(source, bin_size).save(
            output_path / "individual_tiles" / f"{scan.stem}_tile_{tile}.png"
        )


def write_pngs(scan: FileScan, output: pathlib.Path, bin_size: int) -> None:
    """Write the static PNG charts for one input file."""
    # aggregated inputs are bounded by tiles x cycles; Altair's 5000-row cap is too low
    with alt.data_transformers.disable_max_rows():
        plot_n_per_sequence(scan, output / f"{scan.stem}_n_count_distribution.png")
        plot_n_per_cycle_and_tile(
            scan, output / f"{scan.stem}_n_distribution_per_cycle_and_tile.png"
        )
        plot_individual_tiles(scan, output, bin_size)


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
            if scan.reads_with_n == 0:
                logging.warning(f"No N bases found in {scan.path}.")
            kept.append(scan)
    return kept


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running nspector module...")

    bin_size = args.spatial_bin_size
    scans = usable_scans(
        scan_files(
            args.input,
            args.engine,
            args.threads,
            bin_size,
            raw_dir=args.output if getattr(args, "raw", False) else None,
        )
    )
    if not scans:
        logging.warning("No usable input files; nothing to report.")
        return

    for scan in scans:
        logging.debug(
            f"{scan.path.name}: {scan.reads_with_n} of {scan.total_reads} reads have N bases."
        )

    tables = build_tables(scans, bin_size=bin_size)
    tables.write_parquet(args.output)

    if not getattr(args, "no_html", False):
        path = write_report(scans, tables, args.output)
        logging.info(f"Interactive report written to {path}")
    if not getattr(args, "no_multiqc", False):
        write_multiqc(tables, args.output)
    if getattr(args, "png", False):
        for scan in scans:
            if scan.reads_with_n:
                logging.debug(f"Creating PNG charts for file: {scan.path}")
                write_pngs(scan, args.output, bin_size)

    logging.info("Nspector module finished.")
