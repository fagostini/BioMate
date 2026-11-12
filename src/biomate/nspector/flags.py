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

Status is FAIL when the rate is at or above ``fail_rate``, else WARN.
"""

from __future__ import annotations

import polars

from biomate.nspector.tables import NspectorTables

Z_THRESHOLD = 5.0
DARK_TILE_FRACTION = 0.8
# MAD of a normal distribution is 0.6745 sigma; 1.4826 = 1 / 0.6745
MAD_SCALE = 1.4826

FLAG_SCHEMA = {
    "status": polars.Utf8,
    "kind": polars.Utf8,
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
    return (polars.col("rate") - median) / (MAD_SCALE * mad)


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
            tables.tile_cycle.select("read", "lane", "tile", "cycle", "n"),
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
            polars.lit("outlier among the lane's tiles at this cycle").alias("detail"),
        )
    )

    columns = list(FLAG_SCHEMA)
    return (
        polars.concat(
            [f.select(columns).cast(FLAG_SCHEMA) for f in (dark, tiles, cells)]
        )
        .with_columns(polars.col("status").eq("FAIL").alias("_fail"))
        .sort(
            ["_fail", "kind", "read", "lane", "tile", "cycle"],
            descending=[True, False, False, False, False, False],
            nulls_last=True,
        )
        .drop("_fail")
    )


def sample_status(tables: NspectorTables, flags: polars.DataFrame) -> polars.DataFrame:
    """Worst flag status (PASS/WARN/FAIL) per sample over the lanes it ran on."""
    rank = {"PASS": 0, "WARN": 1, "FAIL": 2}
    worst = flags.group_by("read", "lane").agg(
        polars.col("status").replace_strict(rank, return_dtype=polars.Int8).max()
    )
    return (
        tables.sample_lane.select("sample", "read", "lane")
        .join(worst, on=["read", "lane"], how="left")
        .group_by("sample")
        .agg(polars.col("status").fill_null(0).max())
        .with_columns(
            polars.col("status").replace_strict(
                {v: k for k, v in rank.items()}, return_dtype=polars.Utf8
            )
        )
        .sort("sample")
    )
