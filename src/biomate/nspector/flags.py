"""Anomaly flags for nspector, computed from the summary tables (no extra scan).

All thresholds apply to one scale: the fraction of reads with an N at a cycle.
A tile's whole-read rate is its mean per-cycle rate over the cycles that are not
dark, so a tile at 0.01 has, on average, 1% of its reads with an N at each cycle,
and a lane-wide dark cycle is reported once rather than on every tile.

Flag kinds:

* ``dark_cycle``: at least ``DARK_TILE_FRACTION`` of a lane's tiles are at or
  above ``warn_rate`` at that cycle.
* ``tile``: the tile's mean per-cycle rate (dark cycles left out) is at or
  above ``warn_rate``.
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
DARK_TILE_FRACTION = 0.8
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


def _robust_z(group: list[str]) -> polars.Expr:
    """(rate - median) / (1.4826 * MAD) within ``group``; inf when MAD is 0."""
    median = polars.col("rate").median().over(group)
    mad = (polars.col("rate") - median).abs().median().over(group)
    # 0/0 (rate at the median of a zero-spread group) is NaN, which polars
    # orders above every number: a tile at the median is not an outlier
    return ((polars.col("rate") - median) / (MAD_SCALE * mad)).fill_nan(0.0)


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

    tiles = (
        dense.join(dark, on=lane_cycle, how="anti")
        .group_by("read", "lane", "tile")
        .agg(polars.col("rate").mean())
        .with_columns(_robust_z(["read", "lane"]).alias("z"))
        .filter(polars.col("rate") >= warn_rate)
        .with_columns(
            polars.lit("tile").alias("kind"),
            *flowcell,
            polars.lit(None, polars.UInt16).alias("cycle"),
            polars.lit("mean N rate per cycle over the whole read").alias("detail"),
            _status(fail_rate),
        )
    )

    cells = (
        dense.with_columns(_robust_z(lane_cycle).alias("z"))
        .filter((polars.col("z") >= Z_THRESHOLD) & (polars.col("rate") >= warn_rate))
        .join(dark, on=lane_cycle, how="anti")
        .with_columns(_status(fail_rate))
        .join(
            tiles.select("read", "lane", "tile", tile_status="status"),
            on=["read", "lane", "tile"],
            how="left",
        )
        # keep a cell on a flagged tile only when it is worse than the tile flag
        .filter(
            polars.col("tile_status").is_null()
            | ((polars.col("status") == "FAIL") & (polars.col("tile_status") == "WARN"))
        )
        .with_columns(
            polars.lit("tile_cycle").alias("kind"),
            *flowcell,
            polars.lit("outlier among the lane's tiles at this cycle").alias("detail"),
        )
    )

    samples = (
        tables.sample_lane.filter(polars.col("project") != UNDETERMINED)
        .with_columns(polars.col("n_base_rate").alias("rate"))
        .with_columns(_robust_z(["read", "lane"]).alias("z"))
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
            [f.select(columns).cast(FLAG_SCHEMA) for f in (dark, tiles, cells, samples)]
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
    """Fraction of all N bases on a flagged tile, cell or dark cycle.

    High: the Ns come from the flagged events. Low: they are spread over the
    run. ``None`` without N bases. ``tile_cycle.n`` (reads with an N at a
    cycle) counts N bases, so nothing is counted twice.
    """
    cells = tables.tile_cycle.select("read", "lane", "tile", "cycle", "n")
    total = int(cells["n"].sum())
    if not total:
        return None
    keys = {
        "tile": ["read", "lane", "tile"],
        "dark_cycle": ["read", "lane", "cycle"],
        "tile_cycle": ["read", "lane", "tile", "cycle"],
    }
    for kind, on in keys.items():
        hit = flags.filter(polars.col("kind") == kind).select(on).unique()
        cells = cells.join(hit.with_columns(polars.lit(True).alias(kind)), on=on, how="left")
    explained = cells.filter(
        polars.any_horizontal(polars.col(list(keys)).is_not_null())
    )["n"].sum()
    return int(explained) / total
