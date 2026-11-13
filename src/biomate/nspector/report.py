"""Interactive HTML report and MultiQC custom content for nspector."""

from __future__ import annotations

import datetime
import html
import json
import logging
import math
import pathlib
import textwrap
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
# width of the full-row charts, in px (page is 1400 px wide, minus legends)
FULL_WIDTH = 1100
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
# Project / Lane / Read controls: 'All' sums projects or lanes (global view);
# Read = 'All' overlays the reads in line and bar charts and sums them in maps
READ_FILTER = "(read_sel == 'All' || datum.read == read_sel)"
LANE_FILTER = f"(lane_sel == 'All' || datum.lane == lane_sel) && {READ_FILTER}"
SELECTION_FILTER = (
    f"{LANE_FILTER} && (project_sel == 'All' || datum.project == project_sel)"
)


# per read: a base colour and a lighter shade (R1 dark, R2 light, as asked)
READ_SHADES = {
    "N": ["#d62728", "#ff9896", "#8c564b", "#c49c94"],
    "low-Q": ["#4c78a8", "#9ecae9", "#54a24b", "#88d27a"],
}


def _read_color(reads: list[str], title: str | None = "Read") -> alt.Color:
    """Colour by read: R1 in the base blue, R2 in a lighter blue, and so on."""
    return alt.Color("read:N", title=title).scale(
        domain=reads, range=READ_SHADES["low-Q"][: len(reads)]
    )


def _caption(chart: alt.TopLevelMixin, text: str, width: int) -> alt.VConcatChart:
    """Put ``text`` under ``chart`` as a small grey caption, wrapped to ``width`` px."""
    return alt.vconcat(chart).properties(
        title=alt.TitleParams(
            # ~5 px per character at 11 px; Vega titles do not wrap by themselves
            textwrap.wrap(text, max(40, width // 5)),
            orient="bottom",
            anchor="start",
            fontSize=11,
            fontWeight="normal",
            color="#5b6474",
            offset=10,
        )
    )


CAPTIONS = {
    "overview": (
        "Each cell is one sample on one lane, coloured by the metric chosen above "
        "(N rate, N count or low-Q rate) for the selected project and read; with "
        "Read = All the reads are added together. Click a sample to focus the "
        "per-sample charts on it; shift-click adds samples, double-click clears."
    ),
    "sample_cycle": (
        "The chosen metric at each sequencing cycle, one line per sample and read, "
        "added up over the selected lanes; R1 in the base colour, R2 lighter. "
        "Above 200 samples a band shows the median and interquartile range instead. "
        "A spike at one cycle points to a dark cycle. Drag across the chart to pick "
        "a cycle range for the tile map and the spatial map."
    ),
    "n_per_read": (
        "How many N bases each read carries (log scale), per read, for the clicked "
        "samples or all samples in the selection. A few Ns per read are typical "
        "random no-calls; a second hump or a bar at the full read length points to "
        "whole-cycle or whole-read failures."
    ),
    "n_vs_lowq": (
        "Fraction of reads with an N (red) and with a called base below --min-qual "
        "(blue) at each cycle; lighter shades are R2. Ns rising while low-Q stays "
        "flat means the basecaller gave up on those clusters; low-Q rising without "
        "Ns means the signal degraded but bases were still called."
    ),
    "run_length": (
        "For each read with Ns, the length of its longest stretch of consecutive Ns. "
        "Runs of 1 are isolated no-calls; long runs come from bubbles, focus "
        "problems or dark cycles, and a run as long as the read is an all-N read."
    ),
    "first_n": (
        "Where the first N of each read falls: in the first 10 cycles, the last 10 "
        "cycles or in between; all-N reads are counted on their own. Ns early in the "
        "read often point to registration or early chemistry problems, late Ns to "
        "phasing and signal decay."
    ),
    "tile_map": (
        "The flowcell laid out by surface, swath and tile segment, coloured by N "
        "per read in the cycle range picked on the N per cycle chart (all cycles if "
        "none); colours are log-scaled so normal tiles still differ when one is "
        "extreme. Orange outlines mark FAIL tiles (tick Show WARN flags for WARN "
        "too). Click a tile to show it in the spatial map and highlight it in the "
        "tile × cycle heatmap."
    ),
    "spatial": (
        "Where on the tile the N bases sit: N bases counted in square bins of "
        "--spatial-bin-size coordinate units, with totals per column (top) and row "
        "(right), for the clicked tile (or all tiles in the selection) and the "
        "picked cycle range. Ns gathered in one area point to bubbles, debris or a "
        "local focus problem; an even spread points to a run-wide cause."
    ),
    "tile_cycle": (
        "Every tile (rows) at every cycle (columns), coloured by the chosen metric; "
        "empty cells have no N and no low-Q base. A full column is a dark cycle, a "
        "full row a bad tile, and a short streak a local, temporary event such as a "
        "bubble. Colours are log-scaled so normal cells still differ next to an "
        "extreme one; hover for exact values. Orange outlines mark FAIL flags (tick "
        "Show WARN flags for WARN too). Click a row to select its tile."
    ),
    "interop": (
        "One point per tile: its N rate against the InterOp metric picked with the "
        "InterOp control (%PF, cluster density or %≥Q30), coloured by lane. Tiles "
        "with a high N rate and a low %PF or Q30 are likely physical problems on "
        "the flowcell. Click a point to select the tile."
    ),
}


# log-like colours (symlog keeps 0 at the bottom): one extreme tile or cycle no
# longer flattens the rest of the map into one colour; tooltips keep exact values.
# The tiny constant makes rates of 0.1% and 1% clearly different shades.
# unionWith keeps the domain from collapsing to [0, 0] when a selection has no
# N at all; otherwise Vega paints those zeros in the middle colour
ROBUST_SCALE = alt.Scale(
    type="symlog", constant=1e-4, scheme=SCHEME, domain={"unionWith": [0, 1e-6]}
)
# rate metrics as percentages, counts with thousands separators
METRIC_LABELS = (
    "metric == 'count' ? format(datum.value, ',.0f') : format(datum.value, '.2%')"
)
EMPTY_TEXT = "No N bases in this selection"


def _empty_note(
    chart: alt.Chart, field: str, x: float, y: float, text: str = EMPTY_TEXT
) -> alt.Chart:
    """Text shown at (x, y) px when ``sum(field)`` over ``chart``'s data is 0.

    ``chart`` must keep a row for every selected sample/tile even without Ns
    (an aggregate over no rows yields nothing in Vega, so the note would never
    appear on truly empty data).
    """
    return (
        chart.transform_aggregate(total=f"sum({field})")
        .transform_filter("datum.total == 0")
        .mark_text(text=text, fontSize=14, color="#5b6474")
        .encode(x=alt.value(x), y=alt.value(y))
    )


def _log_ticks(maximum: int) -> list[int]:
    """0 and the powers of ten up to ``maximum``, for symlog count axes."""
    ticks, power = [0], 1
    while power <= max(maximum, 1):
        ticks.append(power)
        power *= 10
    return ticks


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
        "project", "read", "lane", "tile", "cycle", "n", "lowq", "total"
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
            .group_by("project", "read", "lane", "tile", "cycle")
            .agg(
                (polars.col("n").sum() / window).alias("n"),
                (polars.col("lowq").sum() / window).alias("lowq"),
                polars.col("total").first(),
            )
            .sort("project", "read", "lane", "tile", "cycle")
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


def with_tile_anchors(
    frame: polars.DataFrame, tables: NspectorTables
) -> polars.DataFrame:
    """Add one zero row (cycle 1) per project/read/lane/tile.

    ``tile_cycle`` is sparse, so a project with no N on a tile at some cycle has
    no row there. When the report sums projects or lanes, every group must still
    contribute its reads to the denominator; the anchors guarantee that.
    """
    anchors = (
        tables.tile_summary.group_by("project", "read", "lane", "tile")
        .agg(polars.col("total").sum())
        .with_columns(
            polars.lit(1, polars.UInt16).alias("cycle"),
            polars.lit(0).alias("n"),
            polars.lit(0).alias("lowq"),
        )
        .select(frame.columns)
    )
    return polars.concat([frame, anchors], how="vertical_relaxed")


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
    projects = sorted(tables.sample_lane["project"].unique().to_list())

    lane_sel = alt.param(
        name="lane_sel",
        value="All",
        bind=alt.binding_select(options=["All", *lanes], name="Lane "),
    )
    project_sel = alt.param(
        name="project_sel",
        value="All",
        bind=alt.binding_select(options=["All", *projects], name="Project "),
    )
    read_options = ["All", *reads] if len(reads) > 1 else reads
    read_sel = alt.param(
        name="read_sel",
        value=read_options[0],
        bind=alt.binding_select(options=read_options, name="Read "),
    )
    read_color = _read_color(reads)
    show_warn = alt.param(
        name="show_warn",
        value=False,
        bind=alt.binding_checkbox(name="Show WARN flags "),
    )
    # FAIL flags are always outlined; WARN only when the checkbox is ticked
    flag_filter = "datum.status == 'FAIL' || show_warn"
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
    pick_sample = alt.selection_point(name="pick_sample", fields=["project", "sample"])
    pick_tile = alt.selection_point(name="pick_tile", fields=["tile"])
    cycle_brush = alt.selection_interval(name="cycle_brush", encodings=["x"])
    # click a read in any read legend to show only it (shift-click adds,
    # double-click shows all); shared by every chart coloured by read
    pick_read = alt.selection_point(name="pick_read", fields=["read"], bind="legend")
    pick_series = alt.selection_point(
        name="pick_series", fields=["series"], bind="legend"
    )

    def shown(selection, opacity: float = 1.0):
        """Opacity ``opacity`` for selected legend entries, 0 (hidden) otherwise."""
        return alt.when(selection).then(alt.value(opacity)).otherwise(alt.value(0))

    color = alt.Color("value:Q", title=_value_title(), scale=ROBUST_SCALE)
    rate_tooltip = [
        alt.Tooltip("n:Q", title="Reads with N", format=","),
        alt.Tooltip("total:Q", title="Reads scanned", format=","),
        alt.Tooltip("rate:Q", title="N rate", format=".3%"),
        alt.Tooltip("lowq_rate:Q", title="Low-Q rate", format=".3%"),
    ]
    rates = {"rate": "datum.n / datum.total", "lowq_rate": "datum.lowq / datum.total"}

    # -- A: sample × lane overview ------------------------------------------
    n_samples = tables.sample_lane.select("project", "sample").unique().height
    label = {"label": "datum.project + ' / ' + datum.sample"}
    overview = (
        alt.Chart(
            tables.sample_lane.select(
                "project", "sample", "read", "lane", "n", "lowq", "total"
            )
        )
        .transform_filter(
            f"{READ_FILTER} && (project_sel == 'All' || datum.project == project_sel)"
        )
        # Read = All: sum the reads of each sample and lane
        .transform_aggregate(
            n="sum(n)",
            lowq="sum(lowq)",
            total="sum(total)",
            groupby=["project", "sample", "lane"],
        )
        .transform_calculate(value=VALUE_EXPR, **rates, **label)
        .mark_rect()
        .encode(
            # at most 8 lanes, so lanes go on Y and samples get the full width
            x=alt.X(
                "label:N",
                title="Project / sample",
                sort="ascending",
                axis=alt.Axis(labelAngle=-45, labelLimit=200),
            ),
            y=alt.Y("lane:O", title="Lane"),
            color=color,
            opacity=alt.when(pick_sample).then(alt.value(1)).otherwise(alt.value(0.3)),
            stroke=alt.when("datum.lane == lane_sel")
            .then(alt.value(FLAG_COLOR))
            .otherwise(alt.value(None)),
            strokeWidth=alt.value(3),
            tooltip=["project:N", "sample:N", "lane:O", *rate_tooltip],
        )
        .add_params(pick_sample)
        .properties(
            title="Samples × lanes (click a sample; outlined = selected lane)",
            width=FULL_WIDTH,
            height=alt.Step(24),
        )
    )

    # -- B: per-sample N rate along the read, with cycle brush --------------
    cycle_base = (
        alt.Chart(
            tables.sample_cycle.select(
                "project", "sample", "read", "lane", "cycle", "n", "lowq", "total"
            )
        )
        .transform_filter(SELECTION_FILTER)
        # Lane = All: one line per sample, summed over its lanes
        .transform_aggregate(
            n="sum(n)",
            lowq="sum(lowq)",
            total="sum(total)",
            groupby=["project", "sample", "read", "cycle"],
        )
        .transform_calculate(
            value=VALUE_EXPR, **rates, label=f"{label['label']} + ' ' + datum.read"
        )
    )
    max_cycle = int(tables.sample_cycle["cycle"].max())
    cycle_ticks = alt.Axis(values=list(range(0, max_cycle + 1, 10)))
    cycle_x = alt.X(
        "cycle:Q", title="Cycle (drag to select a window)", axis=cycle_ticks
    )
    if n_samples <= MANY_SAMPLES:
        sample_cycle = (
            cycle_base.mark_line(interpolate="step-after")
            .encode(
                x=cycle_x,
                y=alt.Y(
                    "value:Q",
                    title=_value_title(),
                    axis=alt.Axis(labelExpr=METRIC_LABELS),
                ),
                detail="label:N",
                color=read_color,
                strokeWidth=alt.when(pick_sample, empty=False)
                .then(alt.value(3))
                .otherwise(alt.value(1.2)),
                opacity=alt.when(~pick_read)
                .then(alt.value(0))
                .when(pick_sample)
                .then(alt.value(0.9))
                .otherwise(alt.value(0.15)),
                tooltip=["project:N", "sample:N", "read:N", "cycle:Q", *rate_tooltip],
            )
            .add_params(cycle_brush, pick_read)
            .properties(
                title="N per cycle, one line per sample and read "
                "(clicked samples drawn thicker)"
            )
        )
    else:
        stats = cycle_base.transform_aggregate(
            median="median(value)",
            q1="q1(value)",
            q3="q3(value)",
            groupby=["read", "cycle"],
        )
        band = (
            stats.mark_area(interpolate="step-after")
            .encode(
                x=cycle_x,
                y=alt.Y("q1:Q", title=_value_title()),
                y2="q3:Q",
                color=read_color,
                opacity=shown(pick_read, 0.3),
                tooltip=[
                    "read:N",
                    "cycle:Q",
                    alt.Tooltip("median:Q", format=".4g"),
                    alt.Tooltip("q1:Q", format=".4g"),
                    alt.Tooltip("q3:Q", format=".4g"),
                ],
            )
            .add_params(cycle_brush, pick_read)
        )
        median = stats.mark_line(interpolate="step-after").encode(
            x="cycle:Q", y="median:Q", color=read_color, opacity=shown(pick_read)
        )
        picked = (
            cycle_base.transform_filter({"param": "pick_sample", "empty": False})
            .mark_line(interpolate="step-after", strokeWidth=3)
            .encode(
                x="cycle:Q",
                y="value:Q",
                detail="label:N",
                color=read_color,
                opacity=shown(pick_read),
                tooltip=["project:N", "sample:N", "read:N", "cycle:Q", *rate_tooltip],
            )
        )
        sample_cycle = alt.layer(band, median, picked).properties(
            title=f"N per cycle: median and IQR of {n_samples} samples per read, "
            "clicked samples drawn thicker"
        )
    sample_cycle = sample_cycle.properties(width=FULL_WIDTH, height=260)

    # -- C: flowcell tile map, restricted to the brushed cycles -------------
    tile_cycle_data, window = window_tile_cycle(tables)
    tile_cycle_data = with_tile_anchors(tile_cycle_data, tables)
    group = ["project", "read", "lane", "tile"]
    tile_base = (
        alt.Chart(tile_cycle_data)
        .transform_filter(SELECTION_FILTER)
        .transform_filter(cycle_brush)
        # per project/read/lane first (total repeats on every cycle row), then per tile
        .transform_aggregate(
            n="sum(n)", lowq="sum(lowq)", total="max(total)", groupby=group
        )
        .transform_aggregate(
            n="sum(n)", lowq="sum(lowq)", total="sum(total)", groupby=["tile"]
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
            )
            .transform_filter(LANE_FILTER)
            .transform_filter(flag_filter)
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
    # size the cells so the map roughly matches the spatial map next to it
    layout = tables.tile_summary.select("surface", "swath", "segment")
    n_segments = max(layout["segment"].n_unique(), 1)
    n_rows = max(layout.select("surface", "swath").unique().height, 1)
    tile_map = alt.layer(tile_map, flagged_tiles).properties(
        title="Flowcell tiles (cycle window from the line chart; click a tile; "
        "orange = flagged)",
        width=alt.Step(tile_step := max(12, min(48, 480 // n_segments))),
        height=alt.Step(max(20, min(60, 480 // n_rows))),
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
        .transform_filter(LANE_FILTER)
        .transform_filter(flag_filter)
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
        .transform_filter(SELECTION_FILTER)
        # each tile's denominator: the reads of every selected project/read/lane on it
        .transform_window(first="row_number()", groupby=group)
        .transform_calculate(group_total="datum.first == 1 ? datum.total : 0")
        .transform_joinaggregate(tile_total="sum(group_total)", groupby=["tile"])
        .transform_aggregate(
            n="sum(n)",
            lowq="sum(lowq)",
            total="max(tile_total)",
            groupby=["tile", "cycle"],
        )
        .transform_filter("datum.n > 0 || datum.lowq > 0")
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
    # clicking a cell selects its tile, like the tile map
    heat_width = max(2, min(8, 1100 * window // max_cycle)) * len(
        range(1, max_cycle + 1, window)
    )
    no_events = _empty_note(
        alt.Chart(tile_cycle_data)
        .transform_filter(SELECTION_FILTER)
        .transform_calculate(events="datum.n + datum.lowq"),
        "events",
        heat_width / 2,
        40,
        "No N bases and no low-Q bases in this selection",
    )
    tile_cycle = alt.layer(
        tile_cycle.add_params(pick_tile).properties(name="tile_cycle_cells"),
        flagged_cells,
        dark_cycles,
        no_events,
    ).properties(
        title="Tile × cycle (empty cells have no N and no low-Q base; "
        "orange = flagged cell or dark cycle)"
        + (f"; cycles merged in windows of {window}" if window > 1 else ""),
        width=alt.Step(cell := max(2, min(8, 1100 * window // max_cycle))),
        # rows at least 6 px tall so a tile can be clicked
        height=alt.Step(max(6, min(14, 700 // max(n_tiles, 1)))),
    )

    # -- E: spatial distribution for the clicked tile -----------------------
    size = tables.bin_size
    spatial = (
        alt.Chart(select_spatial(tables))
        .transform_filter(SELECTION_FILTER)
        .transform_filter(pick_tile)
        # spatial rows are per cycle window; "cycle" is the window's first cycle
        .transform_filter(cycle_brush)
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
        spatial.transform_aggregate(n="sum(n_bases)", groupby=["x_bin", "y_bin"])
        .transform_calculate(**bin_edges)
        .mark_rect()
        .encode(
            x=alt.X("x0:Q", title="X", scale=x_scale),
            x2="x1:Q",
            y=alt.Y("y0:Q", title="Y", scale=y_scale),
            y2="y1:Q",
            color=alt.Color("n:Q", title="N bases", scale=ROBUST_SCALE),
            tooltip=[
                alt.Tooltip("x0:Q", title="X from", format=".0f"),
                alt.Tooltip("y0:Q", title="Y from", format=".0f"),
                alt.Tooltip("n:Q", title="N bases", format=","),
            ],
        )
        .properties(width=480, height=480)
    )
    # anchored tile data keeps a row for every selected tile, so an empty
    # tile/cycle selection still aggregates to 0
    xy = alt.layer(
        xy,
        _empty_note(
            alt.Chart(tile_cycle_data)
            .transform_filter(SELECTION_FILTER)
            .transform_filter(pick_tile)
            .transform_filter(cycle_brush),
            "n",
            240,
            240,
            "No N bases in this tile / cycle selection",
        ),
    )
    top = (
        spatial.transform_aggregate(n="sum(n_bases)", groupby=["x_bin"])
        .transform_calculate(x0=bin_edges["x0"], x1=bin_edges["x1"])
        .mark_rect()
        .encode(
            x=alt.X("x0:Q", title=None, axis=None, scale=x_scale),
            x2="x1:Q",
            y=alt.Y("n:Q", title=None),
        )
        .properties(width=480, height=60)
    )
    right = (
        spatial.transform_aggregate(n="sum(n_bases)", groupby=["y_bin"])
        .transform_calculate(y0=bin_edges["y0"], y1=bin_edges["y1"])
        .mark_rect()
        .encode(
            y=alt.Y("y0:Q", title=None, axis=None, scale=y_scale),
            y2="y1:Q",
            x=alt.X("n:Q", title=None),
        )
        .properties(width=60, height=480)
    )
    spatial_view = alt.vconcat(
        top, alt.hconcat(xy, right, spacing=4), spacing=4
    ).properties(
        title=[
            "Spatial N density: N bases per bin",
            "Tile: clicked on the tile map or tile × cycle (all tiles if none).",
            f"Cycles: windows of {tables.cycle_window} starting in the range dragged "
            "on N per cycle (all if none).",
        ]
    )

    # -- F: N bases per read -------------------------------------------------
    # largest bar either histogram can reach (all samples, lanes and reads)
    n_reads_max = int(
        tables.n_per_read.group_by("read", "n_count")
        .agg(polars.col("n_reads").sum())["n_reads"]
        .max()
        or 1
    )
    # one row per selected sample/read/lane even without Ns: drives the
    # "no N bases" note of the per-read charts
    sample_n = (
        alt.Chart(tables.sample_lane.select("project", "sample", "read", "lane", "n"))
        .transform_filter(SELECTION_FILTER)
        .transform_filter(pick_sample)
    )
    n_per_read = (
        alt.Chart(tables.n_per_read)
        .transform_filter(SELECTION_FILTER)
        .transform_filter(pick_sample)
        .transform_aggregate(n_reads="sum(n_reads)", groupby=["n_count", "read"])
        .mark_bar()
        .encode(
            x=alt.X("n_count:O", title="N bases in read"),
            xOffset="read:N",
            y=alt.Y(
                "n_reads:Q",
                title="Reads",
                axis=alt.Axis(
                    values=_log_ticks(n_reads_max),
                    labelExpr="format(datum.value, ',')",
                ),
            ).scale(type="symlog"),
            color=read_color,
            opacity=shown(pick_read),
            tooltip=["read:N", "n_count:O", alt.Tooltip("n_reads:Q", format=",")],
        )
        .add_params(pick_read)
    )
    n_per_read = alt.layer(
        n_per_read, _empty_note(sample_n, "n", FULL_WIDTH / 2, 120)
    ).properties(
        title="N bases per read (clicked sample, or all)",
        width=FULL_WIDTH,
        height=240,
    )

    # -- G: N vs low-Q per cycle (all samples of the lane) -------------------
    n_series = [f"{kind} {read}" for kind in READ_SHADES for read in reads]
    n_series_colors = [
        READ_SHADES[kind][i % 4] for kind in READ_SHADES for i in range(len(reads))
    ]
    n_vs_lowq = (
        alt.Chart(
            tables.sample_cycle.select(
                "project", "sample", "read", "lane", "cycle", "n", "lowq", "total"
            )
        )
        .transform_filter(SELECTION_FILTER)
        .transform_filter(pick_sample)
        .transform_aggregate(
            n="sum(n)",
            lowq="sum(lowq)",
            total="sum(total)",
            groupby=["cycle", "read"],
        )
        .transform_calculate(
            N="datum.n / datum.total", **{"low-Q": "datum.lowq / datum.total"}
        )
        .transform_fold(["N", "low-Q"], as_=["kind", "fraction"])
        .transform_calculate(series="datum.kind + ' ' + datum.read")
        .mark_line(interpolate="step-after")
        .encode(
            x=alt.X("cycle:Q", title="Cycle", axis=cycle_ticks),
            y=alt.Y("fraction:Q", title="Fraction of reads").axis(format="%"),
            color=alt.Color("series:N", title=None)
            .scale(domain=n_series, range=n_series_colors)
            .legend(orient="top-left"),
            opacity=shown(pick_series),
            tooltip=[
                "cycle:Q",
                "read:N",
                "kind:N",
                alt.Tooltip("fraction:Q", format=".3%"),
            ],
        )
        .add_params(pick_series)
        .properties(
            title="N (no call, red) vs low-Q (called below --min-qual, blue) per "
            "cycle; lighter shade = R2",
            width=FULL_WIDTH,
            height=220,
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
            .transform_filter(READ_FILTER)
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
        .transform_filter(SELECTION_FILTER)
        .transform_filter(pick_sample)
    )
    run_length = (
        shape.transform_aggregate(
            n_reads="sum(n_reads)", groupby=["longest_run", "read"]
        )
        .mark_bar()
        .encode(
            x=alt.X("longest_run:O", title="Longest N run in read"),
            xOffset="read:N",
            y=alt.Y(
                "n_reads:Q",
                title="Reads",
                axis=alt.Axis(
                    values=_log_ticks(n_reads_max),
                    labelExpr="format(datum.value, ',')",
                ),
            ).scale(type="symlog"),
            color=read_color,
            opacity=shown(pick_read),
            tooltip=["read:N", "longest_run:O", alt.Tooltip("n_reads:Q", format=",")],
        )
        .add_params(pick_read)
    )
    run_length = alt.layer(
        run_length, _empty_note(sample_n, "n", FULL_WIDTH / 2, 100)
    ).properties(
        title="Longest N run (1 = isolated Ns; long runs = bubbles, dark cycles)",
        width=FULL_WIDTH,
        height=200,
    )
    first_n = (
        shape.transform_calculate(where="datum.all_n ? 'all-N read' : datum.first_n")
        .transform_aggregate(n_reads="sum(n_reads)", groupby=["where", "read"])
        .mark_bar()
        .encode(
            y=alt.Y(
                "where:N",
                title="First N",
                sort=["first 10 cycles", "middle", "last 10 cycles", "all-N read"],
            ),
            yOffset="read:N",
            x=alt.X("n_reads:Q", title="Reads"),
            color=read_color,
            opacity=shown(pick_read),
            tooltip=["read:N", "where:N", alt.Tooltip("n_reads:Q", format=",")],
        )
        .add_params(pick_read)
    )
    first_n = alt.layer(
        first_n, _empty_note(sample_n, "n", FULL_WIDTH / 2, 70)
    ).properties(title="Where the first N falls", width=FULL_WIDTH, height=140)

    return (
        alt.vconcat(
            _caption(overview, CAPTIONS["overview"], FULL_WIDTH),
            _caption(sample_cycle, CAPTIONS["sample_cycle"], FULL_WIDTH),
            _caption(n_per_read, CAPTIONS["n_per_read"], FULL_WIDTH),
            _caption(n_vs_lowq, CAPTIONS["n_vs_lowq"], FULL_WIDTH),
            _caption(run_length, CAPTIONS["run_length"], FULL_WIDTH),
            _caption(first_n, CAPTIONS["first_n"], FULL_WIDTH),
            # nested concats share colour scales unless told otherwise
            alt.hconcat(
                _caption(
                    tile_map, CAPTIONS["tile_map"], max(tile_step * n_segments, 360)
                ),
                _caption(spatial_view, CAPTIONS["spatial"], 540),
                spacing=48,
            ).resolve_scale(color="independent"),
            _caption(
                tile_cycle,
                CAPTIONS["tile_cycle"],
                cell * len(range(1, max_cycle + 1, window)),
            ),
            *(
                [_caption(interop_view, CAPTIONS["interop"], 420)]
                if interop_view is not None
                else []
            ),
            spacing=36,
        )
        .add_params(project_sel, lane_sel, read_sel, metric, show_warn)
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
            "project": [s.project for s in scans],
            "sample": [s.sample for s in scans],
            "read": [s.read for s in scans],
            "reads scanned": [s.total_reads for s in scans],
            "reads with N": [s.reads_with_n for s in scans],
            "N bases": [sum(t[4] for t in s.tiles) for s in scans],
            "bases": [sum(t[6] for t in s.tiles) for s in scans],
        }
    ).with_columns(
        (polars.col("reads with N") / polars.col("reads scanned")).alias("N rate"),
        (polars.col("N bases") / polars.col("bases")).alias("N per base"),
    )
    worst_tiles = (
        tables.tile_summary.group_by("read", "lane", "tile")
        .agg(
            polars.col("n").sum().alias("reads with N"),
            polars.col("total").sum().alias("reads scanned"),
            polars.col("n_bases").sum().alias("N bases"),
            polars.col("bases").sum(),
        )
        .with_columns(
            (polars.col("reads with N") / polars.col("reads scanned")).alias("N rate"),
            (polars.col("N bases") / polars.col("bases")).alias("N per base"),
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


# number formats of the HTML summary tables
TABLE_FORMATS = {
    "N rate": ".3%",
    "N per base": ".4%",
    "reads with N": ",",
    "reads scanned": ",",
    "N bases": ",",
    "bases": ",",
}

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
  tbody tr:nth-child(even) {{ background: var(--panel); }}
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
      <li>Use the <b>Project</b>, <b>Lane</b>, <b>Read</b> and <b>Metric</b> controls above the charts to switch views. <i>All</i> combines every project or lane (global view); pick <i>Undetermined</i> to see only the reads that were not assigned to a sample. <i>Rate</i> divides by the reads scanned on that tile or sample.</li>
      <li><b>Click a read in a chart legend</b> to show only that read (or series) in the charts that have it; shift-click adds another, double-click shows all again.</li>
      <li><b>Click a sample</b> in the overview to highlight its cycle curve and filter the N-per-read histogram. Shift-click adds more; double-click clears.</li>
      <li><b>Drag on the cycle chart</b> to restrict the flowcell tile map to a cycle window.</li>
      <li><b>Click a tile</b> on the tile map or a row of the tile × cycle heatmap to highlight it in the heatmap and show where on the tile its N reads are.</li>
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
        worst_tiles=_table(worst_tiles, TABLE_FORMATS),
        worst_cycles=_table(worst_cycles, TABLE_FORMATS),
        inputs=_table(inputs, TABLE_FORMATS),
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
