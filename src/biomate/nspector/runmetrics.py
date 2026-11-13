"""Per-tile Illumina InterOp metrics (density, %PF, %>=Q30) for nspector.

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


def _cycle_reads(run_info) -> dict[int, str]:
    """Map every sequencing cycle of a non-index read to R1, R2, ..."""
    reads = run_info.reads()
    cycles: dict[int, str] = {}
    number = 0
    for i in range(reads.size()):
        info = reads[i]
        if info.is_index():
            continue
        number += 1
        for cycle in range(info.first_cycle(), info.last_cycle() + 1):
            cycles[cycle] = f"R{number}"
    return cycles


def read_tile_metrics(run_folder: pathlib.Path) -> polars.DataFrame:
    """Return per read/lane/tile density (K/mm2), %PF and %>=Q30 from InterOp.

    Density and %PF are tile-level and repeated for each read. %>=Q30 sums
    the Q-score histograms over each read's cycles (``QMetricsOut.bin``), or
    uses the Q20/Q30 counts of ``QMetrics2030Out.bin`` when that is what the
    instrument wrote; it is null when neither file exists.
    """
    interop = _import_interop()
    run = interop.read(str(run_folder), valid_to_load=["Tile", "Q", "QCollapsed"])
    if run.tile_metric_set().empty():
        raise RuntimeError(
            f"No tile metrics (InterOp/TileMetricsOut.bin) in {run_folder}."
        )
    cycle_read = _cycle_reads(run.run_info())
    read_names = sorted(set(cycle_read.values()))

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
                )
            )
    # InterOp reports missing values as NaN; store them as nulls
    return polars.DataFrame(rows, schema=INTEROP_SCHEMA, orient="row").fill_nan(None)


def join_tile_interop(
    tables: NspectorTables, metrics: polars.DataFrame
) -> polars.DataFrame:
    """Join InterOp tile metrics to the per read/lane/tile N rate (all samples)."""
    rates = (
        tables.tile_summary.group_by("read", "lane", "tile")
        .agg(polars.col("n", "total").sum())
        .with_columns((polars.col("n") / polars.col("total")).alias("n_rate"))
    )
    return rates.join(
        metrics.cast(
            {k: v for k, v in INTEROP_SCHEMA.items() if k in ("lane", "tile")}
        ),
        on=["read", "lane", "tile"],
        how="left",
    ).sort("read", "lane", "tile")
