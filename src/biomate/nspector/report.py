"""Interactive HTML report and MultiQC custom content for nspector."""

from __future__ import annotations

import datetime
import html
import json
import logging
import math
import pathlib
import warnings

import altair as alt
import polars

from biomate.nspector.flags import FLAG_SCHEMA, sample_status
from biomate.nspector.scan import FileScan
from biomate.nspector.tables import NspectorTables

# Upper bound on spatial rows embedded in the report; beyond this only the
# tiles with the highest N rate keep their spatial map.
MAX_SPATIAL_ROWS = 300_000
# Upper bound on tile x cycle rows embedded in the report; beyond this cycles
# are merged into fixed windows (the Parquet table keeps full detail).
MAX_TILE_CYCLE_ROWS = 200_000
# Above this many samples the cycle chart shows a median/IQR band and only
# draws the clicked samples as lines.
MANY_SAMPLES = 200
# Inlined data above this size triggers a warning (browsers slow down).
SIZE_WARNING_BYTES = 20 * 1024 * 1024

SCHEME = "viridis"
# outline colour for flagged tiles/cells; stands out against both ends of viridis
FLAG_COLOR = "#ff7f0e"
# flags listed in the report header; nspector_flags.parquet keeps all of them
MAX_FLAG_ROWS = 50

# Vega expression: the value plotted for the active metric.
VALUE_EXPR = (
    "metric == 'rate' ? datum.n / datum.total"
    " : metric == 'lowq' ? datum.lowq / datum.total : datum.n"
)
LANE_READ_FILTER = "datum.lane == lane_sel && datum.read == read_sel"


def _value_title() -> str:
    return "N rate / count / low-Q rate"


def select_spatial(tables: NspectorTables, max_rows: int = MAX_SPATIAL_ROWS):
    """Return the spatial rows to embed, trimming to the worst tiles if needed."""
    spatial = tables.spatial
    if spatial.height <= max_rows:
        return spatial
    ranked = (
        tables.tile_summary.group_by("read", "lane", "tile")
        .agg(polars.col("n").sum(), polars.col("total").sum())
        .with_columns((polars.col("n") / polars.col("total")).alias("rate"))
        .join(
            spatial.group_by("read", "lane", "tile").agg(polars.len().alias("rows")),
            on=["read", "lane", "tile"],
        )
        .sort("rate", descending=True)
        .with_columns(polars.col("rows").cum_sum().alias("cum_rows"))
        .filter(polars.col("cum_rows") <= max_rows)
        .select("read", "lane", "tile")
    )
    logging.warning(
        f"Spatial data has {spatial.height} rows; the report keeps the "
        f"{ranked.height} tiles with the highest N rate. Use a larger "
        "--spatial-bin-size to include every tile."
    )
    return spatial.join(ranked, on=["read", "lane", "tile"], how="semi")


def window_tile_cycle(
    tables: NspectorTables, max_rows: int = MAX_TILE_CYCLE_ROWS
) -> tuple[polars.DataFrame, int]:
    """Return the tile x cycle rows to embed and the cycle window size.

    Cycles are merged into windows of ``w`` cycles until the table fits;
    ``n`` becomes the mean per cycle so rates stay comparable.
    """
    frame = tables.tile_cycle.select(
        "read", "lane", "tile", "cycle", "n", "lowq", "total"
    )
    if frame.height <= max_rows:
        return frame, 1
    window = math.ceil(frame.height / max_rows)
    max_cycle = int(frame["cycle"].max())
    while True:
        merged = (
            frame.with_columns(
                ((polars.col("cycle") - 1) // window * window + 1).alias("cycle")
            )
            .group_by("read", "lane", "tile", "cycle")
            .agg(
                (polars.col("n").sum() / window).alias("n"),
                (polars.col("lowq").sum() / window).alias("lowq"),
                polars.col("total").first(),
            )
            .sort("read", "lane", "tile", "cycle")
        )
        # one window per tile is the floor; tiles alone may exceed the budget
        if merged.height <= max_rows or window >= max_cycle:
            break
        window = min(window * 2, max_cycle)
    logging.warning(
        f"Tile x cycle data has {frame.height} rows; the report merges cycles "
        f"into windows of {window}. The Parquet table keeps every cycle."
    )
    return merged, window


def _flags(tables: NspectorTables) -> polars.DataFrame:
    """Return the flags frame, empty when flags were not computed."""
    if tables.flags is None:
        return polars.DataFrame(schema=FLAG_SCHEMA)
    return tables.flags


def _tile_layout_calcs(chart: alt.Chart) -> alt.Chart:
    """Add surface/swath/segment fields derived from the tile number."""
    return chart.transform_calculate(
        surface="datum.tile >= 1000 ? floor(datum.tile / 1000) : 1",
        swath="datum.tile >= 1000 ? floor(datum.tile / 100) % 10 : 1",
        segment="datum.tile >= 1000 ? datum.tile % 100 : datum.tile",
    ).transform_calculate(
        row="'Surface ' + datum.surface + ' / swath ' + datum.swath",
    )


def build_chart(tables: NspectorTables) -> alt.VConcatChart:
    """Build the linked, cross-filtering Vega-Lite report chart."""
    lanes = sorted(tables.sample_lane["lane"].unique().to_list())
    reads = sorted(tables.sample_lane["read"].unique().to_list())

    lane_sel = alt.param(
        name="lane_sel",
        value=lanes[0],
        bind=alt.binding_select(options=lanes, name="Lane "),
    )
    read_sel = alt.param(
        name="read_sel",
        value=reads[0],
        bind=alt.binding_select(options=reads, name="Read "),
    )
    metric = alt.param(
        name="metric",
        value="rate",
        bind=alt.binding_radio(
            options=["rate", "count", "lowq"],
            labels=[
                "N rate (fraction of reads)",
                "N count (reads)",
                "low-Q rate (called bases below --min-qual)",
            ],
            name="Metric ",
        ),
    )
    pick_sample = alt.selection_point(name="pick_sample", fields=["sample"])
    pick_tile = alt.selection_point(name="pick_tile", fields=["tile"])
    cycle_brush = alt.selection_interval(name="cycle_brush", encodings=["x"])

    color = alt.Color("value:Q", title=_value_title()).scale(scheme=SCHEME)
    rate_tooltip = [
        alt.Tooltip("n:Q", title="Reads with N", format=","),
        alt.Tooltip("total:Q", title="Reads scanned", format=","),
        alt.Tooltip("rate:Q", title="N rate", format=".3%"),
        alt.Tooltip("lowq_rate:Q", title="Low-Q rate", format=".3%"),
    ]
    rates = {"rate": "datum.n / datum.total", "lowq_rate": "datum.lowq / datum.total"}

    # -- A: sample × lane overview ------------------------------------------
    n_samples = tables.sample_lane["sample"].n_unique()
    overview = (
        alt.Chart(
            tables.sample_lane.select("sample", "read", "lane", "n", "lowq", "total")
        )
        .transform_filter("datum.read == read_sel")
        .transform_calculate(value=VALUE_EXPR, **rates)
        .mark_rect()
        .encode(
            x=alt.X("lane:O", title="Lane"),
            y=alt.Y("sample:N", title="Sample", sort="ascending"),
            color=color,
            opacity=alt.when(pick_sample).then(alt.value(1)).otherwise(alt.value(0.3)),
            stroke=alt.when("datum.lane == lane_sel")
            .then(alt.value(FLAG_COLOR))
            .otherwise(alt.value(None)),
            strokeWidth=alt.value(3),
            tooltip=["sample:N", "lane:O", *rate_tooltip],
        )
        .add_params(pick_sample)
        .properties(
            title="Samples × lanes (click a sample, outlined = selected lane)",
            width=alt.Step(36),
            height=alt.Step(18) if n_samples < 40 else 600,
        )
    )

    # -- B: per-sample N rate along the read, with cycle brush --------------
    cycle_base = (
        alt.Chart(
            tables.sample_cycle.select(
                "sample", "read", "lane", "cycle", "n", "lowq", "total"
            )
        )
        .transform_filter(LANE_READ_FILTER)
        .transform_calculate(value=VALUE_EXPR, **rates)
    )
    cycle_x = alt.X("cycle:Q", title="Cycle (drag to select a window)")
    if n_samples <= MANY_SAMPLES:
        sample_cycle = (
            cycle_base.mark_line(interpolate="step-after")
            .encode(
                x=cycle_x,
                y=alt.Y("value:Q", title=_value_title()),
                detail="sample:N",
                color=alt.when(pick_sample, empty=False)
                .then(alt.value("#d62728"))
                .otherwise(alt.value("#4c78a8")),
                opacity=alt.when(pick_sample)
                .then(alt.value(0.9))
                .otherwise(alt.value(0.15)),
                tooltip=["sample:N", "cycle:Q", *rate_tooltip],
            )
            .add_params(cycle_brush)
            .properties(title="N per cycle, one line per sample")
        )
    else:
        stats = cycle_base.transform_aggregate(
            median="median(value)",
            q1="q1(value)",
            q3="q3(value)",
            groupby=["cycle"],
        )
        band = (
            stats.mark_area(opacity=0.3, interpolate="step-after", color="#4c78a8")
            .encode(
                x=cycle_x,
                y=alt.Y("q1:Q", title=_value_title()),
                y2="q3:Q",
                tooltip=[
                    "cycle:Q",
                    alt.Tooltip("median:Q", format=".4g"),
                    alt.Tooltip("q1:Q", format=".4g"),
                    alt.Tooltip("q3:Q", format=".4g"),
                ],
            )
            .add_params(cycle_brush)
        )
        median = stats.mark_line(interpolate="step-after", color="#4c78a8").encode(
            x="cycle:Q", y="median:Q"
        )
        picked = (
            cycle_base.transform_filter({"param": "pick_sample", "empty": False})
            .mark_line(interpolate="step-after", color="#d62728")
            .encode(
                x="cycle:Q",
                y="value:Q",
                detail="sample:N",
                tooltip=["sample:N", "cycle:Q", *rate_tooltip],
            )
        )
        sample_cycle = alt.layer(band, median, picked).properties(
            title=f"N per cycle: median and IQR of {n_samples} samples, clicked samples in red"
        )
    sample_cycle = sample_cycle.properties(width=560, height=240)

    # -- C: flowcell tile map, restricted to the brushed cycles -------------
    tile_cycle_data, window = window_tile_cycle(tables)
    tile_base = (
        alt.Chart(tile_cycle_data)
        .transform_filter(LANE_READ_FILTER)
        .transform_filter(cycle_brush)
        .transform_aggregate(
            n="sum(n)", lowq="sum(lowq)", total="max(total)", groupby=["tile"]
        )
        .transform_calculate(value=VALUE_EXPR, **rates)
    )
    flags = _flags(tables)
    flagged_tiles = (
        _tile_layout_calcs(
            alt.Chart(
                flags.filter(polars.col("kind") == "tile").select(
                    "read", "lane", "tile", "status"
                )
            ).transform_filter(LANE_READ_FILTER)
        )
        .mark_rect(filled=False, stroke=FLAG_COLOR, strokeWidth=2)
        .encode(
            x="segment:O",
            y="row:N",
            tooltip=["tile:O", alt.Tooltip("status:N", title="Tile flag")],
        )
    )
    tile_map = (
        _tile_layout_calcs(tile_base)
        .mark_rect()
        .encode(
            x=alt.X("segment:O", title="Tile segment"),
            y=alt.Y("row:N", title=None),
            color=color.title("N per read in selected cycles"),
            opacity=alt.when(pick_tile).then(alt.value(1)).otherwise(alt.value(0.3)),
            tooltip=[
                "tile:O",
                alt.Tooltip("n:Q", title="N bases in window", format=","),
                alt.Tooltip("total:Q", title="Reads scanned", format=","),
                alt.Tooltip("rate:Q", title="N per read", format=".4f"),
                alt.Tooltip("lowq_rate:Q", title="Low-Q per read", format=".4f"),
            ],
        )
        .add_params(pick_tile)
    )
    tile_map = alt.layer(tile_map, flagged_tiles).properties(
        title="Flowcell tiles (cycle window from the line chart; click a tile; "
        "orange = flagged)",
        width=alt.Step(14),
        height=alt.Step(22),
    )

    # -- D: tile × cycle heatmap --------------------------------------------
    n_tiles = tables.tile_cycle["tile"].n_unique()
    max_cycle = int(tables.sample_cycle["cycle"].max())
    cycle_x = alt.X(
        "cycle:O",
        title="Cycle",
        scale=alt.Scale(domain=list(range(1, max_cycle + 1, window))),
        axis=alt.Axis(
            labelExpr=f"(datum.value - 1) % {10 * window} == 0 ? datum.label : ''",
            labelAngle=0,
        ),
    )
    tile_y = alt.Y("tile:O", title="Tile", sort="ascending")
    flag_cells = (
        alt.Chart(
            flags.filter(polars.col("kind").is_in(["tile_cycle", "dark_cycle"])).select(
                "kind", "read", "lane", "tile", "cycle", "status", "rate"
            )
        )
        .transform_filter(LANE_READ_FILTER)
        # map flagged cycles onto the embedded cycle windows
        .transform_calculate(
            cycle=f"floor((datum.cycle - 1) / {window}) * {window} + 1"
        )
    )
    flag_tooltip = [
        "kind:N",
        "status:N",
        "tile:O",
        "cycle:O",
        alt.Tooltip("rate:Q", title="N rate", format=".3%"),
    ]
    flagged_cells = (
        flag_cells.transform_filter("datum.kind == 'tile_cycle'")
        .mark_rect(filled=False, stroke=FLAG_COLOR, strokeWidth=2)
        .encode(x=cycle_x, y=tile_y, tooltip=flag_tooltip)
    )
    dark_cycles = (
        flag_cells.transform_filter("datum.kind == 'dark_cycle'")
        .mark_rect(filled=False, stroke=FLAG_COLOR, strokeWidth=2)
        .encode(x=cycle_x, tooltip=flag_tooltip)
    )
    tile_cycle = (
        alt.Chart(tile_cycle_data)
        .transform_filter(LANE_READ_FILTER)
        .transform_calculate(value=VALUE_EXPR, **rates)
        .mark_rect()
        .encode(
            x=cycle_x,
            y=tile_y,
            color=color,
            opacity=alt.when(pick_tile).then(alt.value(1)).otherwise(alt.value(0.25)),
            tooltip=["tile:O", "cycle:O", *rate_tooltip],
        )
    )
    tile_cycle = alt.layer(tile_cycle, flagged_cells, dark_cycles).properties(
        title="Tile × cycle (empty cells have no N and no low-Q base; "
        "orange = flagged cell or dark cycle)"
        + (f"; cycles merged in windows of {window}" if window > 1 else ""),
        width=alt.Step(max(2, min(8, 1100 * window // max_cycle))),
        height=alt.Step(max(3, min(14, 700 // max(n_tiles, 1)))),
    )

    # -- E: spatial distribution for the clicked tile -----------------------
    size = tables.bin_size
    spatial = (
        alt.Chart(select_spatial(tables))
        .transform_filter(LANE_READ_FILTER)
        .transform_filter(pick_tile)
    )
    bin_edges = {
        "x0": f"datum.x_bin * {size}",
        "x1": f"(datum.x_bin + 1) * {size}",
        "y0": f"datum.y_bin * {size}",
        "y1": f"(datum.y_bin + 1) * {size}",
    }
    # marginals use the same coordinate edges and un-niced data extent as the
    # density map, so their bars line up with its columns and rows
    x_scale = alt.Scale(zero=False, nice=False)
    y_scale = alt.Scale(zero=False, nice=False, reverse=True)
    xy = (
        spatial.transform_aggregate(n="sum(n)", groupby=["x_bin", "y_bin"])
        .transform_calculate(**bin_edges)
        .mark_rect()
        .encode(
            x=alt.X("x0:Q", title="X", scale=x_scale),
            x2="x1:Q",
            y=alt.Y("y0:Q", title="Y", scale=y_scale),
            y2="y1:Q",
            color=alt.Color("n:Q", title="Reads with N").scale(scheme=SCHEME),
            tooltip=[
                alt.Tooltip("x0:Q", title="X from", format=".0f"),
                alt.Tooltip("y0:Q", title="Y from", format=".0f"),
                alt.Tooltip("n:Q", title="Reads with N", format=","),
            ],
        )
        .properties(width=420, height=420)
    )
    top = (
        spatial.transform_aggregate(n="sum(n)", groupby=["x_bin"])
        .transform_calculate(x0=bin_edges["x0"], x1=bin_edges["x1"])
        .mark_rect()
        .encode(
            x=alt.X("x0:Q", title=None, axis=None, scale=x_scale),
            x2="x1:Q",
            y=alt.Y("n:Q", title=None),
        )
        .properties(width=420, height=60)
    )
    right = (
        spatial.transform_aggregate(n="sum(n)", groupby=["y_bin"])
        .transform_calculate(y0=bin_edges["y0"], y1=bin_edges["y1"])
        .mark_rect()
        .encode(
            y=alt.Y("y0:Q", title=None, axis=None, scale=y_scale),
            y2="y1:Q",
            x=alt.X("n:Q", title=None),
        )
        .properties(width=60, height=420)
    )
    spatial_view = alt.vconcat(
        top, alt.hconcat(xy, right, spacing=4), spacing=4
    ).properties(title="Spatial N density (clicked tile, or all tiles of the lane)")

    # -- F: N bases per read -------------------------------------------------
    n_per_read = (
        alt.Chart(tables.n_per_read)
        .transform_filter(LANE_READ_FILTER)
        .transform_filter(pick_sample)
        .transform_aggregate(n_reads="sum(n_reads)", groupby=["n_count"])
        .mark_bar()
        .encode(
            x=alt.X("n_count:O", title="N bases in read"),
            y=alt.Y("n_reads:Q", title="Reads").scale(type="symlog"),
            tooltip=["n_count:O", alt.Tooltip("n_reads:Q", format=",")],
        )
        .properties(
            title="N bases per read (clicked sample, or all)", width=420, height=240
        )
    )

    # -- G: N vs low-Q per cycle (all samples of the lane) -------------------
    n_vs_lowq = (
        alt.Chart(
            tables.sample_cycle.select(
                "sample", "read", "lane", "cycle", "n", "lowq", "total"
            )
        )
        .transform_filter(LANE_READ_FILTER)
        .transform_filter(pick_sample)
        .transform_aggregate(
            n="sum(n)", lowq="sum(lowq)", total="sum(total)", groupby=["cycle"]
        )
        .transform_calculate(
            N="datum.n / datum.total", **{"low-Q": "datum.lowq / datum.total"}
        )
        .transform_fold(["N", "low-Q"], as_=["kind", "fraction"])
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("cycle:Q", title="Cycle"),
            y=alt.Y("fraction:Q", title="Fraction of reads").axis(format="%"),
            color=alt.Color("kind:N", title=None)
            .scale(domain=["N", "low-Q"], range=["#d62728", "#4c78a8"])
            .legend(orient="top-left"),
            tooltip=[
                "cycle:Q",
                "kind:N",
                alt.Tooltip("fraction:Q", format=".3%"),
            ],
        )
        .properties(
            title="N (no call) vs low-Q (called, below --min-qual) per cycle",
            width=420,
            height=200,
        )
    )

    # -- I: N rate against InterOp tile metrics (with --run-folder) ---------
    interop_view = None
    if tables.tile_interop is not None:
        interop_metric = alt.param(
            name="interop_metric",
            value="pct_pf",
            bind=alt.binding_radio(
                options=["pct_pf", "density_k", "pct_q30"],
                labels=["% PF", "Density (K/mm²)", "% ≥ Q30"],
                name="InterOp ",
            ),
        )
        interop_view = (
            alt.Chart(tables.tile_interop)
            .transform_filter("datum.read == read_sel")
            .transform_calculate(
                metric_value="datum[interop_metric]",
                rate="datum.n / datum.total",
            )
            .transform_filter("isValid(datum.metric_value)")
            .mark_circle(size=60)
            .encode(
                x=alt.X("metric_value:Q", title="InterOp metric (see control)").scale(
                    zero=False
                ),
                y=alt.Y("rate:Q", title="N rate (reads with N)").axis(format="%"),
                color=alt.Color("lane:N", title="Lane"),
                opacity=alt.when(pick_tile)
                .then(alt.value(0.9))
                .otherwise(alt.value(0.2)),
                tooltip=[
                    "lane:O",
                    "tile:O",
                    alt.Tooltip("rate:Q", title="N rate", format=".3%"),
                    alt.Tooltip("pct_pf:Q", title="% PF", format=".1f"),
                    alt.Tooltip("density_k:Q", title="Density K/mm²", format=".0f"),
                    alt.Tooltip("pct_q30:Q", title="% ≥ Q30", format=".1f"),
                ],
            )
            .add_params(interop_metric, pick_tile)
            .properties(
                title="N rate vs InterOp tile metric (click a tile)",
                width=420,
                height=240,
            )
        )

    # -- H: N-run shape ------------------------------------------------------
    shape = (
        alt.Chart(tables.n_shape)
        .transform_filter(LANE_READ_FILTER)
        .transform_filter(pick_sample)
    )
    run_length = (
        shape.transform_aggregate(n_reads="sum(n_reads)", groupby=["longest_run"])
        .mark_bar()
        .encode(
            x=alt.X("longest_run:O", title="Longest N run in read"),
            y=alt.Y("n_reads:Q", title="Reads").scale(type="symlog"),
            tooltip=["longest_run:O", alt.Tooltip("n_reads:Q", format=",")],
        )
        .properties(
            title="Longest N run (1 = isolated Ns; long runs = bubbles, dark cycles)",
            width=420,
            height=200,
        )
    )
    first_n = (
        shape.transform_calculate(where="datum.all_n ? 'all-N read' : datum.first_n")
        .transform_aggregate(n_reads="sum(n_reads)", groupby=["where"])
        .mark_bar()
        .encode(
            y=alt.Y(
                "where:N",
                title="First N",
                sort=["first 10 cycles", "middle", "last 10 cycles", "all-N read"],
            ),
            x=alt.X("n_reads:Q", title="Reads"),
            tooltip=["where:N", alt.Tooltip("n_reads:Q", format=",")],
        )
        .properties(title="Where the first N falls", width=200, height=100)
    )

    return (
        alt.vconcat(
            alt.hconcat(overview, sample_cycle),
            alt.hconcat(tile_map, n_per_read),
            alt.hconcat(n_vs_lowq, run_length, first_n),
            tile_cycle,
            *([interop_view] if interop_view is not None else []),
            spatial_view,
            spacing=36,
        )
        .add_params(lane_sel, read_sel, metric)
        .resolve_scale(color="independent")
        .configure_view(fill="#f3f3f3", stroke=None)
        .configure_axis(labelFontSize=11, titleFontSize=12)
    )


def _table(frame: polars.DataFrame, formats: dict[str, str] | None = None) -> str:
    """Render a small polars frame as an HTML table."""
    formats = formats or {}
    head = "".join(f"<th>{html.escape(c)}</th>" for c in frame.columns)
    rows = []
    for row in frame.iter_rows(named=True):
        cells = []
        for col, value in row.items():
            fmt = formats.get(col)
            text = format(value, fmt) if fmt and value is not None else str(value)
            cells.append(f"<td>{html.escape(text)}</td>")
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return (
        f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(rows)}</tbody></table>"
    )


def summary_tables(scans: list[FileScan], tables: NspectorTables, top: int = 10):
    """Return (inputs, worst tiles, worst cycles) frames for the report header."""
    tile_status = (
        _flags(tables)
        .filter(polars.col("kind") == "tile")
        .select("read", "lane", "tile", "status")
    )
    inputs = polars.DataFrame(
        {
            "file": [s.path.name for s in scans],
            "sample": [s.sample for s in scans],
            "read": [s.read for s in scans],
            "reads scanned": [s.total_reads for s in scans],
            "reads with N": [s.reads_with_n for s in scans],
        }
    ).with_columns(
        (polars.col("reads with N") / polars.col("reads scanned")).alias("N rate")
    )
    worst_tiles = (
        tables.tile_summary.group_by("read", "lane", "tile")
        .agg(
            polars.col("n").sum().alias("reads with N"),
            polars.col("total").sum().alias("reads scanned"),
        )
        .with_columns(
            (polars.col("reads with N") / polars.col("reads scanned")).alias("N rate")
        )
        .sort(
            ["N rate", "read", "lane", "tile"], descending=[True, False, False, False]
        )
        .head(top)
        .join(
            tile_status, on=["read", "lane", "tile"], how="left", maintain_order="left"
        )
        .with_columns(polars.col("status").fill_null("PASS"))
    )
    worst_cycles = (
        tables.sample_cycle.group_by("read", "lane", "cycle")
        .agg(
            polars.col("n").sum().alias("reads with N"),
            polars.col("total").sum().alias("reads scanned"),
        )
        .with_columns(
            (polars.col("reads with N") / polars.col("reads scanned")).alias("N rate")
        )
        .sort(
            ["N rate", "read", "lane", "cycle"], descending=[True, False, False, False]
        )
        .head(top)
    )
    return inputs, worst_tiles, worst_cycles


_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>nspector report</title>
<style>
  :root {{ --fg: #1d2330; --muted: #5b6474; --line: #dde1e8; --bg: #ffffff; --panel: #f6f7f9; --accent: #3b5bdb; }}
  body {{ margin: 0; background: var(--bg); color: var(--fg); font: 14px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }}
  header, main {{ max-width: 1400px; margin: 0 auto; padding: 0 16px; }}
  header {{ padding-top: 24px; }}
  h1 {{ font-size: 22px; margin: 0 0 4px; }}
  h2 {{ font-size: 16px; margin: 28px 0 8px; }}
  .meta {{ color: var(--muted); margin: 0; }}
  .stats {{ display: flex; flex-wrap: wrap; gap: 12px; margin: 16px 0; }}
  .stat {{ background: var(--panel); border: 1px solid var(--line); border-radius: 8px; padding: 10px 14px; min-width: 150px; }}
  .stat b {{ display: block; font-size: 20px; font-variant-numeric: tabular-nums; }}
  .stat span {{ color: var(--muted); font-size: 12px; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 16px; }}
  .scroll {{ overflow-x: auto; }}
  table {{ border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }}
  th, td {{ text-align: left; padding: 4px 8px; border-bottom: 1px solid var(--line); white-space: nowrap; }}
  th {{ color: var(--muted); font-weight: 600; font-size: 12px; }}
  .help {{ background: var(--panel); border-left: 3px solid var(--accent); padding: 8px 14px; border-radius: 4px; }}
  .help li {{ margin: 2px 0; }}
  #vis {{ margin: 8px 0 48px; overflow-x: auto; }}
  #controls {{ position: sticky; top: 0; z-index: 2; background: var(--bg); border-bottom: 1px solid var(--line); margin-top: 16px; }}
  #controls .vega-bindings {{ display: flex; flex-wrap: wrap; gap: 4px 24px; padding: 8px 0; }}
  #controls label {{ cursor: pointer; }}
  #controls .vega-bind-name {{ font-weight: 600; }}
</style>
<script type="text/javascript">
{vegaembed}
</script>
</head>
<body>
<header>
  <h1>nspector report</h1>
  <p class="meta">Generated {generated} · {n_files} input file(s) · read length up to {max_len} cycles</p>
  <div class="stats">
    <div class="stat"><b>{total_reads}</b><span>reads scanned</span></div>
    <div class="stat"><b>{n_reads}</b><span>reads with at least one N</span></div>
    <div class="stat"><b>{n_rate}</b><span>overall N rate</span></div>
    <div class="stat"><b>{n_tiles}</b><span>tiles across {n_lanes} lane(s)</span></div>
  </div>
</header>
<main>
  <section><h2>Anomaly flags</h2>{flags}</section>
  <div class="help">
    <ul>
      <li>Use the <b>Lane</b>, <b>Read</b> and <b>Metric</b> controls above the charts to switch views. <i>Rate</i> divides by the reads scanned on that tile or sample.</li>
      <li><b>Click a sample</b> in the overview to highlight its cycle curve and filter the N-per-read histogram. Shift-click adds more; double-click clears.</li>
      <li><b>Drag on the cycle chart</b> to restrict the flowcell tile map to a cycle window.</li>
      <li><b>Click a tile</b> on the tile map to highlight it in the heatmap and show where on the tile its N reads are.</li>
      <li><i>Low-Q rate</i> counts reads with a called base (not N) below <code>--min-qual</code>. High N with normal low-Q means the basecaller gave up; high low-Q means degraded but called.</li>
    </ul>
  </div>
  <div id="controls"></div>
  <div id="vis"></div>
  <div class="grid">
    <section><h2>Tiles with the highest N rate</h2><div class="scroll">{worst_tiles}</div></section>
    <section><h2>Cycles with the highest N rate</h2><div class="scroll">{worst_cycles}</div></section>
  </div>
  <section><h2>Inputs</h2><div class="scroll">{inputs}</div></section>
  <p class="meta" style="margin:24px 0">Full tables: <code>nspector_*.parquet</code> in this directory.</p>
</main>
<script type="text/javascript">
  const spec = {spec};
  vegaEmbed('#vis', spec, {{"mode": "vega-lite", "renderer": "canvas", "actions": {{"export": true, "source": false, "compiled": false, "editor": false}}}})
    .then((result) => {{
      window.nspectorView = result.view;
      const bindings = document.querySelector('#vis .vega-bindings');
      if (bindings) document.getElementById('controls').appendChild(bindings);
    }})
    .catch(console.error);
</script>
</body>
</html>
"""


def _vegaembed_bundle() -> str:
    """Return the vega + vega-lite + vega-embed JS bundle (offline)."""
    from altair.utils._importers import import_vl_convert, vl_version_for_vl_convert

    vlc = import_vl_convert()
    return vlc.javascript_bundle(vl_version=vl_version_for_vl_convert())


def write_report(
    scans: list[FileScan], tables: NspectorTables, output: pathlib.Path
) -> pathlib.Path:
    """Write the self-contained interactive report to ``output/nspector_report.html``."""
    # row counts are bounded by MAX_SPATIAL_ROWS / MAX_TILE_CYCLE_ROWS instead
    with alt.data_transformers.disable_max_rows():
        # pick_tile is shared on purpose by the tile map and the InterOp scatter
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message="Automatically deduplicated selection parameter"
            )
            spec = json.dumps(build_chart(tables).to_dict(), separators=(",", ":"))
    if len(spec) > SIZE_WARNING_BYTES:
        logging.warning(
            f"The report embeds {len(spec) / 2**20:.0f} MB of data and may be slow "
            "to open. Consider a larger --spatial-bin-size."
        )
    # keep "</script>" in string data from closing the script element
    spec = spec.replace("</", "<\\/")

    inputs, worst_tiles, worst_cycles = summary_tables(scans, tables)
    flags = _flags(tables)
    if flags.height:
        n_fail = int((flags["status"] == "FAIL").sum())
        more = (
            f"; showing the first {MAX_FLAG_ROWS}, all in nspector_flags.parquet"
            if flags.height > MAX_FLAG_ROWS
            else ""
        )
        flags_html = (
            f'<p class="meta">{flags.height} flag(s): {n_fail} FAIL, '
            f"{flags.height - n_fail} WARN{more}.</p>"
            '<div class="scroll">'
            + _table(flags.head(MAX_FLAG_ROWS), {"rate": ".3%", "z": ".1f"})
            + "</div>"
        )
    else:
        flags_html = "<p class='meta'>No tile, cell or dark-cycle flags.</p>"
    total_reads = int(tables.sample_lane["total"].sum())
    n_reads = int(tables.sample_lane["n"].sum())
    page = _PAGE.format(
        vegaembed=_vegaembed_bundle(),
        flags=flags_html,
        spec=spec,
        generated=datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
        n_files=len(scans),
        max_len=int(tables.sample_cycle["cycle"].max()),
        total_reads=f"{total_reads:,}",
        n_reads=f"{n_reads:,}",
        n_rate=f"{n_reads / total_reads:.3%}" if total_reads else "n/a",
        n_tiles=tables.tile_summary.select("lane", "tile").unique().height,
        n_lanes=tables.sample_lane["lane"].n_unique(),
        worst_tiles=_table(
            worst_tiles, {"N rate": ".3%", "reads with N": ",", "reads scanned": ","}
        ),
        worst_cycles=_table(
            worst_cycles, {"N rate": ".3%", "reads with N": ",", "reads scanned": ","}
        ),
        inputs=_table(
            inputs, {"N rate": ".3%", "reads with N": ",", "reads scanned": ","}
        ),
    )
    path = output / "nspector_report.html"
    path.write_text(page, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# MultiQC custom content
# ---------------------------------------------------------------------------


def _row_label(sample: str, lane: int, read: str) -> str:
    return f"{sample} L{lane:03d} {read}"


def write_multiqc(tables: NspectorTables, output: pathlib.Path) -> list[pathlib.Path]:
    """Write MultiQC custom-content files under ``output/multiqc``.

    Run ``multiqc <output>/multiqc`` (or include the folder in a larger
    MultiQC run) to pick them up.
    """
    folder = output / "multiqc"
    folder.mkdir(parents=True, exist_ok=True)
    written = []

    # tile heatmap: rows = sample/lane/read, columns = tiles, % reads with N
    tile = tables.tile_summary.with_columns(
        polars.struct("sample", "lane", "read")
        .map_elements(
            lambda r: _row_label(r["sample"], r["lane"], r["read"]),
            return_dtype=polars.Utf8,
        )
        .alias("row")
    )
    ycats = sorted(tile["row"].unique().to_list())
    xcats = [str(t) for t in sorted(tile["tile"].unique().to_list())]
    lookup = {
        (r["row"], str(r["tile"])): round(100 * r["n_rate"], 4)
        for r in tile.select("row", "tile", "n_rate").iter_rows(named=True)
    }
    heatmap = {
        "id": "nspector_tile_rate",
        "section_name": "nspector: N rate per tile",
        "description": "Percentage of reads with at least one N base on each flowcell tile.",
        "plot_type": "heatmap",
        "pconfig": {
            "id": "nspector_tile_rate_heatmap",
            "title": "nspector: % reads with N per tile",
            "xlab": "Tile",
            "ylab": "Sample / lane / read",
        },
        "xcats": xcats,
        "ycats": ycats,
        "data": [[lookup.get((y, x), 0.0) for x in xcats] for y in ycats],
    }
    path = folder / "nspector_tile_rate_mqc.json"
    path.write_text(json.dumps(heatmap, indent=1))
    written.append(path)

    # cycle line graph: one series per sample/lane/read
    series: dict[str, dict[int, float]] = {}
    for r in tables.sample_cycle.iter_rows(named=True):
        label = _row_label(r["sample"], r["lane"], r["read"])
        series.setdefault(label, {})[int(r["cycle"])] = round(100 * r["n_rate"], 4)
    linegraph = {
        "id": "nspector_cycle_rate",
        "section_name": "nspector: N rate per cycle",
        "description": "Percentage of reads with an N base at each sequencing cycle.",
        "plot_type": "linegraph",
        "pconfig": {
            "id": "nspector_cycle_rate_plot",
            "title": "nspector: % reads with N per cycle",
            "xlab": "Cycle",
            "ylab": "% reads with N",
            "ymin": 0,
        },
        "data": series,
    }
    path = folder / "nspector_cycle_rate_mqc.json"
    path.write_text(json.dumps(linegraph, indent=1))
    written.append(path)

    flags = _flags(tables)
    if flags.height:
        table = {
            "id": "nspector_flags",
            "section_name": "nspector: anomaly flags",
            "description": "Tiles, tile-cycle cells and dark cycles with high N rates "
            "(see nspector_flags.parquet).",
            "plot_type": "table",
            "pconfig": {"id": "nspector_flags_table", "title": "nspector: flags"},
            "data": {
                f"{i + 1:04d} {r['read']} L{r['lane']:03d}": {
                    k: v
                    for k, v in r.items()
                    # JSON has no inf (z is inf when the lane's MAD is 0)
                    if v is not None and not (isinstance(v, float) and math.isinf(v))
                }
                for i, r in enumerate(flags.iter_rows(named=True))
            },
        }
        path = folder / "nspector_flags_mqc.json"
        path.write_text(json.dumps(table, indent=1))
        written.append(path)

    # general statistics: one row per sample (all lanes and reads)
    per_sample = (
        tables.sample_lane.group_by("sample")
        .agg(
            polars.col("n").sum(),
            polars.col("n_bases").sum(),
            polars.col("lowq").sum(),
            polars.col("total").sum(),
        )
        .with_columns(
            (100 * polars.col("n") / polars.col("total")).alias("pct"),
            (100 * polars.col("lowq") / polars.col("total")).alias("pct_lowq"),
        )
        .join(sample_status(tables, flags), on="sample", how="left")
        .sort("sample")
    )
    lines = [
        "# id: 'nspector'",
        "# section_name: 'nspector'",
        "# plot_type: 'generalstats'",
        "# pconfig:",
        "#     - pct_reads_with_n:",
        "#         title: '% N reads'",
        "#         description: 'Percentage of reads with at least one N base (nspector)'",
        "#         min: 0",
        "#         max: 100",
        "#         suffix: '%'",
        "#         format: '{:,.3f}'",
        "#     - n_bases:",
        "#         title: 'N bases'",
        "#         description: 'Total N base calls (nspector)'",
        "#         format: '{:,.0f}'",
        "#     - pct_reads_lowq:",
        "#         title: '% low-Q reads'",
        "#         description: 'Percentage of reads with a called base below --min-qual (nspector)'",
        "#         min: 0",
        "#         max: 100",
        "#         suffix: '%'",
        "#         format: '{:,.3f}'",
        "#     - n_status:",
        "#         title: 'N status'",
        "#         description: 'Worst nspector flag on the sample lanes (PASS/WARN/FAIL)'",
        "Sample\tpct_reads_with_n\tn_bases\tpct_reads_lowq\tn_status",
    ]
    lines += [
        f"{r['sample']}\t{r['pct']:.4f}\t{r['n_bases']}\t{r['pct_lowq']:.4f}\t{r['status']}"
        for r in per_sample.iter_rows(named=True)
    ]
    path = folder / "nspector_general_stats_mqc.tsv"
    path.write_text("\n".join(lines) + "\n")
    written.append(path)
    return written
