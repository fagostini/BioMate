"""Tidy summary tables built from nspector scans.

Every table carries both the number of N-containing reads (``n``) and the
number of reads scanned (``total``) so that rates can be derived at any
aggregation level without losing the denominator.
"""

from __future__ import annotations

import dataclasses
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
    bin_size: int

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


def _frame(scan: FileScan, rows: list[tuple], schema: dict) -> polars.DataFrame:
    """Build a frame from a FileScan counter list, tagged with sample and read."""
    return polars.DataFrame(rows, schema=schema, orient="row").with_columns(
        polars.lit(scan.sample).alias("sample"),
        polars.lit(scan.read).alias("read"),
    )


LANE_TILE = {"lane": polars.UInt8, "tile": polars.UInt32}


def _concat(scans: list[FileScan], field: str, schema: dict) -> polars.DataFrame:
    return polars.concat([_frame(s, getattr(s, field), schema) for s in scans])


def _rate(frame: polars.DataFrame) -> polars.DataFrame:
    return frame.with_columns((polars.col("n") / polars.col("total")).alias("n_rate"))


def build_tables(scans: list[FileScan], bin_size: int = 1000) -> NspectorTables:
    """Combine per-file scan counters into the tidy nspector tables.

    * ``sample_lane``: reads with N per sample/read/lane.
    * ``sample_cycle``: reads with an N at each cycle, per sample/read/lane
      (dense over 1..read length).
    * ``tile_summary``: per sample/read/lane/tile counts, N bases and rate.
    * ``tile_cycle``: per read/lane/tile/cycle, summed over samples (the
      flowcell is physical, so tile defects are shared by every sample);
      sparse, only cycles with at least one N.
    * ``n_per_read``: histogram of N bases per read.
    * ``spatial``: N-containing reads per ``bin_size`` x ``bin_size``
      coordinate bin, per read/lane/tile, summed over samples; sparse.
    """
    tiles = _concat(
        scans,
        "tiles",
        {
            **LANE_TILE,
            "total": polars.UInt64,
            "n": polars.UInt64,
            "n_bases": polars.UInt64,
        },
    )
    cycles = _concat(
        scans,
        "tile_cycle",
        {**LANE_TILE, "cycle": polars.UInt16, "n": polars.UInt64},
    )
    max_len = polars.DataFrame(
        {
            "sample": [s.sample for s in scans],
            "read": [s.read for s in scans],
            "max_len": [s.max_len for s in scans],
        },
        schema={"sample": polars.Utf8, "read": polars.Utf8, "max_len": polars.UInt16},
    )

    # -- sample level -----------------------------------------------------
    sample_lane = (
        _rate(tiles.group_by(KEYS).agg(polars.col("n", "n_bases", "total").sum()))
        .select(*KEYS, "n", "n_bases", "total", "n_rate")
        .sort(KEYS)
    )

    cycle_grid = (
        sample_lane.select(*KEYS, "total")
        .join(
            max_len.group_by("sample", "read").agg(polars.col("max_len").max()),
            on=["sample", "read"],
        )
        .with_columns(
            polars.int_ranges(1, polars.col("max_len").cast(polars.Int64) + 1)
            .cast(polars.List(polars.UInt16))
            .alias("cycle")
        )
        .explode("cycle")
    )
    sample_cycle = (
        _rate(
            cycle_grid.join(
                cycles.group_by(*KEYS, "cycle").agg(polars.col("n").sum()),
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
            tiles.group_by(*KEYS, "tile").agg(polars.col("n", "n_bases", "total").sum())
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

    flowcell_totals = tiles.group_by("read", "lane", "tile").agg(
        polars.col("total").sum()
    )
    tile_cycle = _rate(
        cycles.group_by("read", "lane", "tile", "cycle")
        .agg(polars.col("n").sum())
        .join(flowcell_totals, on=["read", "lane", "tile"], how="left")
    ).sort(["read", "lane", "tile", "cycle"])

    # -- read level -------------------------------------------------------
    n_per_read = (
        _concat(
            scans,
            "n_per_read",
            {"lane": polars.UInt8, "n_count": polars.UInt16, "n_reads": polars.UInt64},
        )
        .group_by(*KEYS, "n_count")
        .agg(polars.col("n_reads").sum())
        .sort([*KEYS, "n_count"])
    )

    spatial = (
        _concat(
            scans,
            "spatial",
            {
                **LANE_TILE,
                "x_bin": polars.UInt32,
                "y_bin": polars.UInt32,
                "n": polars.UInt64,
            },
        )
        .group_by("read", "lane", "tile", "x_bin", "y_bin")
        .agg(polars.col("n").sum())
        .sort(["read", "lane", "tile", "x_bin", "y_bin"])
    )

    return NspectorTables(
        sample_lane=sample_lane,
        sample_cycle=sample_cycle,
        tile_summary=tile_summary,
        tile_cycle=tile_cycle,
        n_per_read=n_per_read,
        spatial=spatial,
        bin_size=bin_size,
    )
