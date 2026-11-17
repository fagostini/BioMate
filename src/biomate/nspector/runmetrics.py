"""Illumina InterOp metrics for nspector: per tile (density, %PF, %>=Q30,
% occupied, phasing) and per cycle (intensity, focus, error rate, no-calls).

Needs the optional ``interop`` package (``pip install 'biomate[interop]'``),
imported lazily so the rest of nspector works without it.
"""

from __future__ import annotations

import pathlib

import polars

from biomate.nspector.tables import NspectorTables

INTEROP_SCHEMA = {
    "read": polars.Utf8,
    "lane": polars.UInt8,
    "tile": polars.UInt32,
    "density_k": polars.Float64,
    "pct_pf": polars.Float64,
    "pct_q30": polars.Float64,
    "pct_occupied": polars.Float64,
    "pct_phasing": polars.Float64,
    "pct_prephasing": polars.Float64,
}

# per read/lane/tile/cycle; the report averages the tiles of a lane
CYCLE_SCHEMA = {
    "read": polars.Utf8,
    "lane": polars.UInt8,
    "tile": polars.UInt32,
    "cycle": polars.UInt16,
    "intensity": polars.Float64,
    "fwhm": polars.Float64,
    "error_rate": polars.Float64,
    "pct_no_call": polars.Float64,
}


def _import_interop():
    try:
        import interop
    except ImportError as err:
        raise RuntimeError(
            "--run-folder needs the Illumina interop package: "
            "pip install 'biomate[interop]'"
        ) from err
    return interop


def _read_names(run_info) -> dict[int, str]:
    """Map the run's read numbers (index reads left out) to R1, R2, ..."""
    reads = run_info.reads()
    names: dict[int, str] = {}
    for i in range(reads.size()):
        if not reads[i].is_index():
            names[reads[i].number()] = f"R{len(names) + 1}"
    return names


def _cycle_reads(run_info) -> dict[int, str]:
    """Map every sequencing cycle of a non-index read to R1, R2, ..."""
    reads = run_info.reads()
    names = _read_names(run_info)
    cycles: dict[int, str] = {}
    for i in range(reads.size()):
        info = reads[i]
        if info.is_index():
            continue
        for cycle in range(info.first_cycle(), info.last_cycle() + 1):
            cycles[cycle] = names[info.number()]
    return cycles


def read_tile_metrics(run_folder: pathlib.Path) -> polars.DataFrame:
    """Return per read/lane/tile density (K/mm2), %PF, %>=Q30, % occupied and
    phasing from InterOp.

    Density and %PF are tile-level and repeated for each read. %>=Q30 sums
    the Q-score histograms over each read's cycles (``QMetricsOut.bin``), or
    uses the Q20/Q30 counts of ``QMetrics2030Out.bin`` when that is what the
    instrument wrote; it is null when neither file exists. % occupied
    (patterned flowcells, ``ExtendedTileMetricsOut.bin``) and the phasing and
    prephasing rates are null when the instrument did not write them.
    """
    interop = _import_interop()
    run = interop.read(
        str(run_folder), valid_to_load=["Tile", "ExtendedTile", "Q", "QCollapsed"]
    )
    if run.tile_metric_set().empty():
        raise RuntimeError(
            f"No tile metrics (InterOp/TileMetricsOut.bin) in {run_folder}."
        )
    cycle_read = _cycle_reads(run.run_info())
    number_read = _read_names(run.run_info())
    read_names = sorted(number_read.values())
    extended = run.extended_tile_metric_set()
    occupied = {
        (e.lane(), e.tile()): e.percent_occupied()
        for e in (extended.at(i) for i in range(extended.size()))
    }

    q30: dict[tuple[str, int, int], list[int]] = {}

    def add(lane: int, tile: int, cycle: int, over: int, total: int) -> None:
        read = cycle_read.get(cycle)
        if read is not None:
            acc = q30.setdefault((read, lane, tile), [0, 0])
            acc[0] += over
            acc[1] += total

    q_set = run.q_metric_set()
    if not q_set.empty():
        index = q_set.index_for_q_value(30)
        for i in range(q_set.size()):
            q = q_set.at(i)
            add(
                q.lane(),
                q.tile(),
                q.cycle(),
                q.total_over_qscore(index),
                q.sum_qscore(),
            )
    else:
        collapsed = run.q_collapsed_metric_set()
        for i in range(collapsed.size()):
            q = collapsed.at(i)
            add(q.lane(), q.tile(), q.cycle(), q.q30(), q.total())

    rows = []
    tiles = run.tile_metric_set()
    for i in range(tiles.size()):
        t = tiles.at(i)
        phasing = {
            number_read.get(r.read()): (r.percent_phasing(), r.percent_prephasing())
            for r in t.read_metrics()
        }
        for read in read_names:
            over, total = q30.get((read, t.lane(), t.tile()), (0, 0))
            rows.append(
                (
                    read,
                    t.lane(),
                    t.tile(),
                    t.cluster_density() / 1000,
                    t.percent_pf(),
                    100 * over / total if total else None,
                    occupied.get((t.lane(), t.tile())),
                    *phasing.get(read, (None, None)),
                )
            )
    # InterOp reports missing values as NaN; store them as nulls
    return polars.DataFrame(rows, schema=INTEROP_SCHEMA, orient="row").fill_nan(None)


def cycle_frame(run) -> polars.DataFrame:
    """Per read/lane/tile/cycle InterOp metrics from loaded ``run`` metrics.

    ``intensity`` and ``fwhm`` (spot width, larger = out of focus) average
    the channels; ``cycle`` counts within the read. Index reads are left out.
    """
    interop = _import_interop()
    table = interop.imaging(run)
    names = table.dtype.names or ()
    if not len(table):
        return polars.DataFrame(schema=CYCLE_SCHEMA)
    frame = polars.DataFrame({name: table[name] for name in names}).fill_nan(None)

    def mean_of(prefix: str) -> polars.Expr:
        columns = [name for name in names if name.startswith(prefix)]
        if not columns:
            return polars.lit(None, polars.Float64)
        return polars.mean_horizontal(columns)

    return (
        frame.select(
            polars.col("Read")
            .cast(polars.Int64)
            .replace_strict(
                _read_names(run.run_info()), default=None, return_dtype=polars.Utf8
            )
            .alias("read"),
            polars.col("Lane").alias("lane"),
            polars.col("Tile").alias("tile"),
            polars.col("Cycle Within Read").alias("cycle"),
            mean_of("P90/").alias("intensity"),
            mean_of("Fwhm/").alias("fwhm"),
            mean_of("Error Rate").alias("error_rate"),
            mean_of("% No Calls").alias("pct_no_call"),
        )
        .drop_nulls("read")
        # one row per tile and cycle already; the mean only merges duplicates
        .group_by("read", "lane", "tile", "cycle")
        .mean()
        .cast(CYCLE_SCHEMA)
        .select(list(CYCLE_SCHEMA))
        .sort("read", "lane", "tile", "cycle")
    )


def read_cycle_metrics(run_folder: pathlib.Path) -> polars.DataFrame:
    """Return :func:`cycle_frame` for the InterOp files of ``run_folder``."""
    interop = _import_interop()
    return cycle_frame(
        interop.read(str(run_folder), valid_to_load=interop.load_imaging_metrics())
    )


def join_tile_interop(
    tables: NspectorTables, metrics: polars.DataFrame
) -> polars.DataFrame:
    """Join InterOp tile metrics to the per read/lane/tile N rates (all samples)."""
    rates = (
        tables.tile_summary.group_by("read", "lane", "tile")
        .agg(polars.col("n", "total", "n_bases", "bases").sum())
        .with_columns(
            (polars.col("n") / polars.col("total")).alias("n_rate"),
            (polars.col("n_bases") / polars.col("bases")).alias("n_base_rate"),
        )
    )
    return rates.join(
        metrics.cast(
            {k: v for k, v in INTEROP_SCHEMA.items() if k in ("lane", "tile")}
        ),
        on=["read", "lane", "tile"],
        how="left",
    ).sort("read", "lane", "tile")
