"""Tidy summary tables built from nspector scans.

Every table carries both the number of N-containing reads (``n``) and the
number of reads scanned (``total``) so that rates can be derived at any
aggregation level without losing the denominator.
"""

from __future__ import annotations

import dataclasses
import logging
import pathlib

import polars

from biomate.nspector.scan import FileScan

KEYS = ["sample", "read", "lane"]


@dataclasses.dataclass
class NspectorTables:
    """All summary tables produced by an nspector run."""

    sample_lane: polars.DataFrame
    sample_cycle: polars.DataFrame
    tile_summary: polars.DataFrame
    tile_cycle: polars.DataFrame
    n_per_read: polars.DataFrame
    spatial: polars.DataFrame
    spatial_grid: dict[str, float]

    def items(self) -> list[tuple[str, polars.DataFrame]]:
        """Return (name, frame) pairs for every table."""
        return [
            (f.name, getattr(self, f.name))
            for f in dataclasses.fields(self)
            if isinstance(getattr(self, f.name), polars.DataFrame)
        ]

    def write_parquet(self, output: pathlib.Path) -> list[pathlib.Path]:
        """Write every table as ``nspector_<name>.parquet`` and return the paths."""
        paths = []
        for name, frame in self.items():
            path = output / f"nspector_{name}.parquet"
            frame.write_parquet(path)
            paths.append(path)
        return paths


def tile_layout(tile: polars.Expr) -> list[polars.Expr]:
    """Split Illumina tile numbers into surface, swath and segment columns.

    4+ digit tiles follow ``<surface><swath><segment>`` (e.g. 1101, 2278).
    Shorter tile numbers are treated as surface 1, swath 1.
    """
    text = tile.cast(polars.Utf8)
    long = text.str.len_chars() >= 4
    return [
        polars.when(long)
        .then(text.str.slice(0, 1).cast(polars.Int16))
        .otherwise(1)
        .alias("surface"),
        polars.when(long)
        .then(text.str.slice(1, 1).cast(polars.Int16))
        .otherwise(1)
        .alias("swath"),
        polars.when(long)
        .then(text.str.slice(2).cast(polars.Int16))
        .otherwise(tile.cast(polars.Int16))
        .alias("segment"),
    ]


def reads_frame(scan: FileScan) -> polars.DataFrame:
    """One row per N-containing read with its coordinates and N cycles."""
    return polars.DataFrame(
        {
            "lane": scan.lanes,
            "tile": scan.tiles,
            "x": scan.xs,
            "y": scan.ys,
            "cycles": scan.cycles,
        },
        schema={
            "lane": polars.UInt8,
            "tile": polars.UInt16,
            "x": polars.UInt32,
            "y": polars.UInt32,
            "cycles": polars.List(polars.UInt16),
        },
    ).with_columns(
        polars.lit(scan.sample).alias("sample"),
        polars.lit(scan.read).alias("read"),
        polars.col("cycles").list.len().cast(polars.UInt16).alias("n_count"),
    )


def totals_frame(scan: FileScan) -> polars.DataFrame:
    """One row per (lane, tile) with the number of reads scanned."""
    return polars.DataFrame(
        {
            "lane": [lane for lane, _ in scan.tile_totals],
            "tile": [tile for _, tile in scan.tile_totals],
            "total": list(scan.tile_totals.values()),
        },
        schema={"lane": polars.UInt8, "tile": polars.UInt16, "total": polars.UInt64},
    ).with_columns(
        polars.lit(scan.sample).alias("sample"),
        polars.lit(scan.read).alias("read"),
        polars.lit(scan.max_len).cast(polars.UInt16).alias("max_len"),
    )


def _rate(frame: polars.DataFrame) -> polars.DataFrame:
    return frame.with_columns((polars.col("n") / polars.col("total")).alias("n_rate"))


def build_tables(scans: list[FileScan], spatial_bins: int = 50) -> NspectorTables:
    """Combine per-file scans into the tidy nspector tables.

    * ``sample_lane``: reads with N per sample/read/lane.
    * ``sample_cycle``: reads with an N at each cycle, per sample/read/lane
      (dense over 1..read length).
    * ``tile_summary``: per sample/read/lane/tile counts, N bases and rate.
    * ``tile_cycle``: per read/lane/tile/cycle, summed over samples (the
      flowcell is physical, so tile defects are shared by every sample);
      sparse, only cycles with at least one N.
    * ``n_per_read``: histogram of N bases per read.
    * ``spatial``: N-containing reads binned on a ``spatial_bins`` grid per
      read/lane/tile, summed over samples; sparse.
    """
    reads = polars.concat([reads_frame(s) for s in scans])
    totals = polars.concat([totals_frame(s) for s in scans])
    exploded = reads.explode("cycles").rename({"cycles": "cycle"})

    # -- sample level -----------------------------------------------------
    sample_totals = totals.group_by(KEYS).agg(
        polars.col("total").sum(), polars.col("max_len").max()
    )
    sample_lane = (
        _rate(
            sample_totals.join(
                reads.group_by(KEYS).agg(
                    polars.len().cast(polars.UInt64).alias("n"),
                    polars.col("n_count").sum().cast(polars.UInt64).alias("n_bases"),
                ),
                on=KEYS,
                how="left",
            ).with_columns(polars.col("n", "n_bases").fill_null(0))
        )
        .select(*KEYS, "n", "n_bases", "total", "n_rate")
        .sort(KEYS)
    )

    cycle_grid = sample_totals.with_columns(
        polars.int_ranges(1, polars.col("max_len").cast(polars.Int64) + 1)
        .cast(polars.List(polars.UInt16))
        .alias("cycle")
    ).explode("cycle")
    sample_cycle = (
        _rate(
            cycle_grid.join(
                exploded.group_by(*KEYS, "cycle").agg(
                    polars.len().cast(polars.UInt64).alias("n")
                ),
                on=[*KEYS, "cycle"],
                how="left",
            ).with_columns(polars.col("n").fill_null(0))
        )
        .select(*KEYS, "cycle", "n", "total", "n_rate")
        .sort([*KEYS, "cycle"])
    )

    # -- tile level -------------------------------------------------------
    tile_summary = (
        _rate(
            totals.join(
                reads.group_by(*KEYS, "tile").agg(
                    polars.len().cast(polars.UInt64).alias("n"),
                    polars.col("n_count").sum().cast(polars.UInt64).alias("n_bases"),
                ),
                on=[*KEYS, "tile"],
                how="left",
            ).with_columns(polars.col("n", "n_bases").fill_null(0))
        )
        .with_columns(
            polars.when(polars.col("n") > 0)
            .then(polars.col("n_bases") / polars.col("n"))
            .otherwise(0.0)
            .alias("mean_n_per_affected_read"),
            *tile_layout(polars.col("tile")),
        )
        .select(
            *KEYS,
            "tile",
            "surface",
            "swath",
            "segment",
            "n",
            "n_bases",
            "total",
            "n_rate",
            "mean_n_per_affected_read",
        )
        .sort([*KEYS, "tile"])
    )

    flowcell_totals = totals.group_by("read", "lane", "tile").agg(
        polars.col("total").sum()
    )
    tile_cycle = _rate(
        exploded.group_by("read", "lane", "tile", "cycle")
        .agg(polars.len().cast(polars.UInt64).alias("n"))
        .join(flowcell_totals, on=["read", "lane", "tile"], how="left")
    ).sort(["read", "lane", "tile", "cycle"])

    # -- read level -------------------------------------------------------
    n_per_read = (
        reads.group_by(*KEYS, "n_count")
        .agg(polars.len().cast(polars.UInt64).alias("n_reads"))
        .sort([*KEYS, "n_count"])
    )

    spatial, grid = _spatial(reads, spatial_bins)

    return NspectorTables(
        sample_lane=sample_lane,
        sample_cycle=sample_cycle,
        tile_summary=tile_summary,
        tile_cycle=tile_cycle,
        n_per_read=n_per_read,
        spatial=spatial,
        spatial_grid=grid,
    )


def _spatial(
    reads: polars.DataFrame, bins: int
) -> tuple[polars.DataFrame, dict[str, float]]:
    """Bin N-containing reads on a ``bins`` × ``bins`` grid per tile."""
    if bins < 1:
        raise ValueError("spatial bins must be >= 1")
    if reads.is_empty():
        empty = polars.DataFrame(
            schema={
                "read": polars.Utf8,
                "lane": polars.UInt8,
                "tile": polars.UInt16,
                "x_bin": polars.UInt16,
                "y_bin": polars.UInt16,
                "n": polars.UInt64,
            }
        )
        return empty, {
            "x_min": 0.0,
            "y_min": 0.0,
            "x_width": 1.0,
            "y_width": 1.0,
            "bins": bins,
        }
    x_min, x_max = int(reads["x"].min()), int(reads["x"].max())
    y_min, y_max = int(reads["y"].min()), int(reads["y"].max())
    x_width = max((x_max - x_min + 1) / bins, 1.0)
    y_width = max((y_max - y_min + 1) / bins, 1.0)
    grid = {
        "x_min": float(x_min),
        "y_min": float(y_min),
        "x_width": x_width,
        "y_width": y_width,
        "bins": bins,
    }
    logging.debug(f"Spatial grid: {grid}")
    spatial = (
        reads.with_columns(
            ((polars.col("x") - x_min) / x_width)
            .floor()
            .clip(0, bins - 1)
            .cast(polars.UInt16)
            .alias("x_bin"),
            ((polars.col("y") - y_min) / y_width)
            .floor()
            .clip(0, bins - 1)
            .cast(polars.UInt16)
            .alias("y_bin"),
        )
        .group_by("read", "lane", "tile", "x_bin", "y_bin")
        .agg(polars.len().cast(polars.UInt64).alias("n"))
        .sort(["read", "lane", "tile", "x_bin", "y_bin"])
    )
    return spatial, grid
