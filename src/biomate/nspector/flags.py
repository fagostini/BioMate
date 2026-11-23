"""Anomaly flags for nspector, computed from the summary tables (no extra scan).

All thresholds apply to one scale: the fraction of reads with an N at a cycle.
A tile's whole-read rate is its mean per-cycle rate over the cycles that are not
dark, so a tile at 0.01 has, on average, 1% of its reads with an N at each cycle,
and a lane-wide dark cycle is reported once rather than on every tile.

Flag kinds:

* ``dark_cycle``: at least ``DARK_TILE_FRACTION`` of a lane's tiles are at or
  above ``warn_rate`` at that cycle. Not raised when ``WHOLE_READ_FRACTION`` of
  the read's cycles qualify: that is bad tiles or a bad lane.
* ``tile``: the tile's mean per-cycle rate (dark cycles left out) is at or
  above ``warn_rate``.
* ``lane``: at least ``LANE_TILE_FRACTION`` of a lane's tiles (two or more)
  qualify for a tile flag. Reported once; tiles and cells on that lane are
  listed only when they are worse than the lane flag.
* ``tile_cycle``: one tile at one cycle is a robust outlier among the lane's
  tiles at that cycle (z >= ``Z_THRESHOLD``) and at or above ``warn_rate``.
  Cells already explained by a dark cycle, or by a tile flag at least as
  severe as the cell, are not repeated.
* ``sample``: the sample's N per sequenced base (its mean per-cycle rate) is a
  robust outlier among the lane's samples (z >= ``Z_THRESHOLD``) and at or
  above ``warn_rate``. Undetermined reads are left out. The only kind that
  names a ``project`` and ``sample``; the others are flowcell-level.

Status is FAIL when the rate is at or above ``fail_rate``, else WARN.
"""

from __future__ import annotations

import polars

from biomate.nspector.scan import UNDETERMINED
from biomate.nspector.tables import KEYS, NspectorTables

Z_THRESHOLD = 5.0
# half the tiles: above it the median tile is dark, so the cells of that cycle
# are no longer outliers and only a dark-cycle flag can report them
DARK_TILE_FRACTION = 0.5
LANE_TILE_FRACTION = 0.8
# a read with this share of its cycles dark has bad tiles or a bad lane, not
# cycle events
WHOLE_READ_FRACTION = 0.8
# MAD of a normal distribution is 0.6745 sigma; 1.4826 = 1 / 0.6745
MAD_SCALE = 1.4826

FLAG_SCHEMA = {
    "status": polars.Utf8,
    "kind": polars.Utf8,
    "project": polars.Utf8,
    "sample": polars.Utf8,
    "read": polars.Utf8,
    "lane": polars.UInt8,
    "tile": polars.UInt32,
    "cycle": polars.UInt16,
    "rate": polars.Float64,
    "z": polars.Float64,
    "detail": polars.Utf8,
}


def _status(fail_rate: float) -> polars.Expr:
    return (
        polars.when(polars.col("rate") >= fail_rate)
        .then(polars.lit("FAIL"))
        .otherwise(polars.lit("WARN"))
        .alias("status")
    )


def _robust_z(group: list[str], warn_rate: float) -> polars.Expr:
    """(rate - median) / max(1.4826 * MAD, warn_rate / Z_THRESHOLD) within ``group``.

    The floor on the scale means a z of ``Z_THRESHOLD`` needs a rate at least
    ``warn_rate`` above the group's median. Without it a group with almost no
    spread (a few near-equal samples, or tiles all at 0) scores any small
    difference as a huge outlier, or as inf.
    """
    median = polars.col("rate").median().over(group)
    mad = (polars.col("rate") - median).abs().median().over(group)
    scale = polars.max_horizontal(MAD_SCALE * mad, polars.lit(warn_rate / Z_THRESHOLD))
    return (polars.col("rate") - median) / scale


def _not_covered(
    frame: polars.DataFrame, cover: polars.DataFrame, on: list[str]
) -> polars.DataFrame:
    """Drop the ``frame`` flags explained by a ``cover`` flag at least as severe."""
    return (
        frame.join(cover.select(*on, _cover="status"), on=on, how="left")
        .filter(
            polars.col("_cover").is_null()
            | ((polars.col("status") == "FAIL") & (polars.col("_cover") == "WARN"))
        )
        .drop("_cover")
    )


def dense_tile_cycle(tables: NspectorTables) -> polars.DataFrame:
    """Per read/lane/tile/cycle N rate, including cycles with no N (rate 0)."""
    cycles = tables.sample_cycle.group_by("read").agg(
        polars.col("cycle").max().alias("max_cycle")
    )
    totals = tables.tile_summary.group_by("read", "lane", "tile").agg(
        polars.col("total").sum()
    )
    grid = (
        totals.join(cycles, on="read")
        .with_columns(
            polars.int_ranges(1, polars.col("max_cycle").cast(polars.Int64) + 1)
            .cast(polars.List(polars.UInt16))
            .alias("cycle")
        )
        .explode("cycle")
        .drop("max_cycle")
    )
    return (
        grid.join(
            # flags are flowcell-level: sum the projects sharing a tile
            tables.tile_cycle.group_by("read", "lane", "tile", "cycle").agg(
                polars.col("n").sum()
            ),
            on=["read", "lane", "tile", "cycle"],
            how="left",
        )
        .with_columns(polars.col("n").fill_null(0))
        .with_columns((polars.col("n") / polars.col("total")).alias("rate"))
    )


def build_flags(
    tables: NspectorTables, warn_rate: float = 0.01, fail_rate: float = 0.05
) -> polars.DataFrame:
    """Return every flag as one frame (schema ``FLAG_SCHEMA``), FAIL first."""
    dense = dense_tile_cycle(tables)
    lane_cycle = ["read", "lane", "cycle"]

    flowcell = [
        polars.lit(None, polars.Utf8).alias("project"),
        polars.lit(None, polars.Utf8).alias("sample"),
    ]

    dark = (
        dense.group_by(lane_cycle)
        .agg(
            (polars.col("rate") >= warn_rate).mean().alias("frac"),
            (polars.col("n").sum() / polars.col("total").sum()).alias("rate"),
            polars.len().alias("tiles"),
        )
        .filter(polars.col("frac") >= DARK_TILE_FRACTION)
        .with_columns(
            polars.lit("dark_cycle").alias("kind"),
            *flowcell,
            polars.lit(None, polars.UInt32).alias("tile"),
            polars.lit(None, polars.Float64).alias("z"),
            polars.format(
                "{}% of {} tiles at or above the warn rate",
                (100 * polars.col("frac")).round(0).cast(polars.Int64),
                polars.col("tiles"),
            ).alias("detail"),
            _status(fail_rate),
        )
    )

    # dark on most of the read: leave it to the tile and lane flags, which then
    # count every cycle, instead of listing nearly all cycles as dark
    whole = (
        dark.group_by("read", "lane")
        .agg(polars.len().alias("dark"))
        .join(
            dense.group_by("read", "lane").agg(polars.col("cycle").n_unique()),
            on=["read", "lane"],
        )
        .filter(polars.col("dark") >= WHOLE_READ_FRACTION * polars.col("cycle"))
    )
    dark = dark.join(whole, on=["read", "lane"], how="anti")

    tile_rates = (
        dense.join(dark, on=lane_cycle, how="anti")
        .group_by("read", "lane", "tile")
        .agg(polars.col("rate").mean())
        .with_columns(_robust_z(["read", "lane"], warn_rate).alias("z"))
    )
    lanes = (
        tile_rates.group_by("read", "lane")
        .agg(
            (polars.col("rate") >= warn_rate).mean().alias("frac"),
            polars.col("rate").mean(),
            polars.len().alias("tiles"),
        )
        .filter((polars.col("frac") >= LANE_TILE_FRACTION) & (polars.col("tiles") > 1))
        .with_columns(
            polars.lit("lane").alias("kind"),
            *flowcell,
            polars.lit(None, polars.UInt32).alias("tile"),
            polars.lit(None, polars.UInt16).alias("cycle"),
            polars.lit(None, polars.Float64).alias("z"),
            polars.format(
                "{}% of {} tiles at or above the warn rate",
                (100 * polars.col("frac")).round(0).cast(polars.Int64),
                polars.col("tiles"),
            ).alias("detail"),
            _status(fail_rate),
        )
    )
    all_tiles = tile_rates.filter(polars.col("rate") >= warn_rate).with_columns(
        polars.lit("tile").alias("kind"),
        *flowcell,
        polars.lit(None, polars.UInt16).alias("cycle"),
        polars.lit("mean N rate per cycle over the whole read").alias("detail"),
        _status(fail_rate),
    )
    tiles = _not_covered(all_tiles, lanes, ["read", "lane"])

    cells = dense.with_columns(_robust_z(lane_cycle, warn_rate).alias("z")).filter(
        (polars.col("z") >= Z_THRESHOLD) & (polars.col("rate") >= warn_rate)
    )
    cells = cells.join(dark, on=lane_cycle, how="anti").with_columns(_status(fail_rate))
    # keep a cell only when it is worse than the flag on its tile or lane
    cells = _not_covered(cells, all_tiles, ["read", "lane", "tile"])
    cells = _not_covered(cells, lanes, ["read", "lane"]).with_columns(
        polars.lit("tile_cycle").alias("kind"),
        *flowcell,
        polars.lit("outlier among the lane's tiles at this cycle").alias("detail"),
    )

    samples = (
        tables.sample_lane.filter(polars.col("project") != UNDETERMINED)
        .with_columns(polars.col("n_base_rate").alias("rate"))
        .with_columns(_robust_z(["read", "lane"], warn_rate).alias("z"))
        .filter((polars.col("z") >= Z_THRESHOLD) & (polars.col("rate") >= warn_rate))
        .with_columns(
            polars.lit("sample").alias("kind"),
            polars.lit(None, polars.UInt32).alias("tile"),
            polars.lit(None, polars.UInt16).alias("cycle"),
            polars.lit("outlier among the lane's samples (N per base)").alias("detail"),
            _status(fail_rate),
        )
    )

    columns = list(FLAG_SCHEMA)
    return (
        polars.concat(
            [
                f.select(columns).cast(FLAG_SCHEMA)
                for f in (dark, lanes, tiles, cells, samples)
            ]
        )
        .with_columns(polars.col("status").eq("FAIL").alias("_fail"))
        .sort(
            ["_fail", "kind", "read", "lane", "tile", "cycle"],
            descending=[True, False, False, False, False, False],
            nulls_last=True,
        )
        .drop("_fail")
    )


RANK = {"PASS": 0, "WARN": 1, "FAIL": 2}


def flag_rank(tables: NspectorTables, flags: polars.DataFrame) -> polars.DataFrame:
    """``sample_lane`` with ``status``: the worst flag rank (see ``RANK``) per row.

    Flowcell-level flags apply to every sample on their read/lane; ``sample``
    flags only to the sample they name.
    """
    rank = polars.col("status").replace_strict(RANK, return_dtype=polars.Int8).max()
    own = polars.col("kind") == "sample"
    return (
        tables.sample_lane.join(
            flags.filter(~own).group_by("read", "lane").agg(rank.alias("lane_rank")),
            on=["read", "lane"],
            how="left",
        )
        .join(
            flags.filter(own).group_by(KEYS).agg(rank.alias("own_rank")),
            on=KEYS,
            how="left",
        )
        .with_columns(
            polars.max_horizontal(
                polars.col("lane_rank", "own_rank").fill_null(0)
            ).alias("status")
        )
    )


def status_name() -> polars.Expr:
    """Turn a ``status`` rank column back into PASS/WARN/FAIL."""
    return polars.col("status").replace_strict(
        {v: k for k, v in RANK.items()}, return_dtype=polars.Utf8
    )


def sample_status(tables: NspectorTables, flags: polars.DataFrame) -> polars.DataFrame:
    """Worst flag status (PASS/WARN/FAIL) per sample over the lanes it ran on."""
    return (
        flag_rank(tables, flags)
        .group_by("sample")
        .agg(polars.col("status").max())
        .with_columns(status_name())
        .sort("sample")
    )


def flagged_share(tables: NspectorTables, flags: polars.DataFrame) -> float | None:
    """Fraction of all N bases on a flagged lane, tile, cell or dark cycle.

    High: the Ns come from the flagged events. Low: they are spread over the
    run. ``None`` without N bases. ``tile_cycle.n`` (reads with an N at a
    cycle) counts N bases, so nothing is counted twice.
    """
    cells = tables.tile_cycle.select("read", "lane", "tile", "cycle", "n")
    total = int(cells["n"].sum())
    if not total:
        return None
    keys = {
        "lane": ["read", "lane"],
        "tile": ["read", "lane", "tile"],
        "dark_cycle": ["read", "lane", "cycle"],
        "tile_cycle": ["read", "lane", "tile", "cycle"],
    }
    marks = [f"_{kind}" for kind in keys]
    for (kind, on), mark in zip(keys.items(), marks):
        hit = flags.filter(polars.col("kind") == kind).select(on).unique()
        cells = cells.join(
            hit.with_columns(polars.lit(True).alias(mark)), on=on, how="left"
        )
    explained = cells.filter(polars.any_horizontal(polars.col(marks).is_not_null()))[
        "n"
    ].sum()
    return int(explained) / total


# tiles named in one row of the collapsed flag table before "and N more"
MAX_TILES_NAMED = 6


def collapse_flags(flags: polars.DataFrame) -> polars.DataFrame:
    """Flags for display: one row per run of cycles, tiles side by side.

    A run of consecutive cycles becomes one row (``cycles`` is ``"31"`` or
    ``"45–52"``), and the ``tile_cycle`` flags of a lane with the same cycles
    become one row naming their tiles (``tile`` is text). ``events`` counts
    the flags in the row; ``rate`` and ``z`` are their largest and ``status``
    their worst. FAIL first, like ``build_flags``.
    """
    keys = ["kind", "project", "sample", "read", "lane"]
    first, last = polars.col("_from"), polars.col("_to")
    tiles = polars.col("tile").drop_nulls().unique().sort().cast(polars.Utf8)
    runs = (
        flags.sort(*keys, "tile", "cycle")
        # consecutive cycles share cycle - position
        .with_columns(
            (
                polars.col("cycle").cast(polars.Int64)
                - polars.int_range(polars.len()).over(*keys, "tile")
            ).alias("_run")
        )
        .group_by(*keys, "tile", "_run")
        .agg(
            polars.col("status").eq("FAIL").any().alias("_fail"),
            polars.col("cycle").min().alias("_from"),
            polars.col("cycle").max().alias("_to"),
            polars.col("rate").max(),
            polars.col("z").max(),
            polars.col("detail").first(),
            polars.len().alias("events"),
        )
        # only tile x cycle cells are merged across tiles
        .with_columns(
            polars.when(polars.col("kind") == "tile_cycle")
            .then(None)
            .otherwise("tile")
            .alias("_own_tile")
        )
    )
    return (
        runs.group_by(*keys, "_from", "_to", "_own_tile")
        .agg(
            polars.col("_fail").any(),
            tiles.head(MAX_TILES_NAMED).str.join(", ").alias("tile"),
            (tiles.len().cast(polars.Int64) - MAX_TILES_NAMED).alias("_more"),
            polars.col("rate").max(),
            polars.col("z").max(),
            polars.col("detail").first(),
            polars.col("events").sum(),
        )
        .with_columns(
            polars.when("_fail")
            .then(polars.lit("FAIL"))
            .otherwise(polars.lit("WARN"))
            .alias("status"),
            polars.when(first == last)
            .then(first.cast(polars.Utf8))
            .otherwise(polars.format("{}–{}", first, last))
            .alias("cycles"),
            polars.when(polars.col("_more") > 0)
            .then(polars.format("{} and {} more", "tile", "_more"))
            # no tile at all (dark cycle, lane, sample) is null, not ""
            .otherwise(polars.col("tile").replace("", None))
            .alias("tile"),
        )
        .sort(
            ["_fail", "kind", "read", "lane", "_from", "tile"],
            descending=[True, False, False, False, False, False],
            nulls_last=True,
        )
        .select("status", *keys, "tile", "cycles", "events", "rate", "z", "detail")
    )


def removal_costs(tables: NspectorTables, flags: polars.DataFrame) -> list[dict]:
    """What removing the flagged tiles and trailing dark cycles would cost.

    One dict per action with ``action``, ``what``, ``lost`` (fraction of the
    reads or bases given up) and ``n_removed`` (fraction of all N bases that
    go with them). Excluding a tile drops its clusters in every read. Only
    dark cycles that run to the end of a read can be trimmed: demultiplexers
    cannot leave out cycles in the middle of a read.
    """
    costs = []
    total_n = int(tables.sample_lane["n_bases"].sum())
    if not total_n:
        return costs

    bad = flags.filter(polars.col("kind") == "tile").select("lane", "tile").unique()
    if bad.height:
        tiles = tables.tile_summary.join(bad, on=["lane", "tile"], how="semi")
        names = (
            bad.sort("lane", "tile")
            .group_by("lane", maintain_order=True)
            .agg(polars.col("tile").cast(polars.Utf8).str.join(", "))
        )
        costs.append(
            {
                "action": f"Exclude the {bad.height} flagged tile(s)",
                "what": "; ".join(f"lane {lane}: {t}" for lane, t in names.iter_rows()),
                "lost": tiles["total"].sum() / tables.tile_summary["total"].sum(),
                "n_removed": tiles["n_bases"].sum() / total_n,
            }
        )

    dark = flags.filter(polars.col("kind") == "dark_cycle")
    lengths = tables.sample_cycle.group_by("read").agg(polars.col("cycle").max())
    for read, length in lengths.sort("read").iter_rows():
        cycles = set(dark.filter(polars.col("read") == read)["cycle"].to_list())
        first = length + 1
        while first - 1 in cycles:
            first -= 1
        if first > length:
            continue
        tail = tables.sample_cycle.filter(
            (polars.col("read") == read) & (polars.col("cycle") >= first)
        )
        bases = tables.sample_lane.filter(polars.col("read") == read)["bases"].sum()
        costs.append(
            {
                "action": f"Trim {read} after cycle {first - 1}",
                "what": f"dark cycles {first}–{length} run to the end of the read",
                # sample_cycle.n counts the reads with an N at the cycle: N bases
                "lost": (tail["total"].sum() - tail["n"].sum()) / bases,
                "n_removed": tail["n"].sum() / total_n,
            }
        )
    return costs


def _some(items: list[str], limit: int = 5) -> str:
    """Join ``items`` with commas, naming at most ``limit`` of them."""
    more = len(items) - limit
    return ", ".join(items[:limit]) + (f" and {more} more" if more > 0 else "")


def verdict(flags: polars.DataFrame, share: float | None) -> str:
    """Sum up the collapsed flags (``collapse_flags``) in a few plain sentences."""
    if not flags.height:
        return "No anomalies: no lane, tile, cycle or sample was flagged."

    def of(kind: str) -> list[dict]:
        return flags.filter(polars.col("kind") == kind).to_dicts()

    def where(row: dict) -> str:
        return f"lane {row['lane']} {row['read']}"

    def fails(rows: list[dict]) -> str:
        count = sum(r["status"] == "FAIL" for r in rows)
        return f" ({count} FAIL)" if count else ""

    parts = []
    if rows := of("lane"):
        names = _some([f"{where(r)} ({r['status']})" for r in rows])
        parts.append(f"Whole lane elevated: {names}.")
    dark: dict[str, list[str]] = {}
    for row in of("dark_cycle"):
        dark.setdefault(where(row), []).append(row["cycles"])
    parts += [f"Dark cycles on {k}: {_some(v)}." for k, v in dark.items()]
    if rows := of("tile"):
        worst = max(rows, key=lambda r: r["rate"])
        parts.append(
            f"{len(rows)} flagged tile(s){fails(rows)}; the worst is "
            f"{worst['tile']} on {where(worst)}, at {worst['rate']:.2%} per cycle."
        )
    if rows := of("tile_cycle"):
        events = sum(r["events"] for r in rows)
        worst = max(r["rate"] for r in rows)
        parts.append(
            f"{events} tile × cycle event(s) at {len(rows)} cycle position(s), "
            f"up to {worst:.0%} of a tile's reads."
        )
    if rows := of("sample"):
        names = dict.fromkeys(f"{r['project']} / {r['sample']}" for r in rows)
        parts.append(f"Outlier sample(s): {_some(list(names))}.")
    if share is not None:
        parts.append(
            f"Flagged lanes, tiles and cycles hold {share:.0%} of all N bases."
        )
    return " ".join(parts)


def read_quality(tables: NspectorTables) -> polars.DataFrame:
    """Per read: share of reads with a poly-G tail, mean Phred of the called
    bases and mean Phred of the called bases next to an N (null without Ns).
    """
    total = polars.col("total").sum()
    return (
        tables.sample_lane.group_by("read")
        .agg(
            (polars.col("polyg").sum() / total).alias("polyg_rate"),
            (
                polars.col("qual_sum").sum()
                / (polars.col("bases").sum() - polars.col("n_bases").sum())
            ).alias("mean_q"),
            (
                polars.col("flank_qual_sum").sum() / polars.col("flank_bases").sum()
            ).alias("flank_q"),
        )
        .fill_nan(None)
        .sort("read")
    )


def summary(
    tables: NspectorTables, flags: polars.DataFrame, params: dict | None = None
) -> dict:
    """The run in one JSON-ready dict: status, verdict, counts, rates and costs.

    ``status`` is the worst flag (PASS without flags). Written as
    ``nspector_summary.json`` for pipelines and LIMS.
    """
    lanes = tables.sample_lane
    counts = dict(flags["status"].value_counts().iter_rows())
    share = flagged_share(tables, flags)
    rates = [
        (polars.col("n").sum() / polars.col("total").sum()).alias("n_rate"),
        (polars.col("n_bases").sum() / polars.col("bases").sum()).alias("n_per_base"),
    ]
    return {
        "status": "FAIL"
        if counts.get("FAIL")
        else "WARN"
        if counts.get("WARN")
        else "PASS",
        "verdict": verdict(collapse_flags(flags), share),
        "flags": {status: counts.get(status, 0) for status in ("FAIL", "WARN")},
        "flags_by_kind": dict(flags["kind"].value_counts().sort("kind").iter_rows()),
        "flagged_n_share": share,
        "reads": int(lanes["total"].sum()),
        "reads_with_n": int(lanes["n"].sum()),
        **lanes.select(rates).row(0, named=True),
        "polyg_rate": int(lanes["polyg"].sum()) / max(int(lanes["total"].sum()), 1),
        "by_read": {
            row.pop("read"): row
            for row in lanes.group_by("read")
            .agg(rates)
            .join(read_quality(tables), on="read")
            .sort("read")
            .to_dicts()
        },
        "removal_costs": removal_costs(tables, flags),
        "parameters": params or {},
    }
