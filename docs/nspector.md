# Nspector

This module inspects one or more FASTQ files for N-base content and shows how those bases are distributed across lanes, samples, flowcell tiles and sequencing cycles.
It helps diagnose quality issues tied to specific physical locations on the flowcell, to particular cycles, or to individual samples.

All inputs are combined into a single run-level result:

- **an interactive HTML report**: one self-contained file that works offline;
- **tidy Parquet tables** for your own analysis;
- **MultiQC custom content**, so the results can appear in a facility MultiQC report.

Static PNG charts are still available with `--png`.

## Usage

The minimal command is `biomate nspector --input <FASTQ_FILE> [<FASTQ_FILE> ...]`.

```bash
biomate nspector \
    --input Run1/SampleA_S1_L00{1,2}_R{1,2}_001.fastq.gz Run1/SampleB_S2_L00{1,2}_R{1,2}_001.fastq.gz \
    --output ./nspector_results --engine rust
```

## Options

- `--input`: one or more FASTQ files (`.fastq`, `.fq`, optionally gzipped). Required.
- `--output`: output directory, created if missing (default: `./output`).
- `--engine`: `python` (default) or `rust`. `rust` uses the native `biomate_kernels` scanner and falls back to Python if it is not installed. Both engines produce identical results.
- `--threads`: extra decompression threads (Python engine) or parallel files (Rust engine).
- `--spatial-bins`: grid size per axis for the per-tile spatial maps (default: `50`). Lower it if the report gets large.
- `--png`: also write the static PNG charts for each input file.
- `--no-html`: skip the HTML report.
- `--no-multiqc`: skip the MultiQC files.

## Counts and rates

Nspector reads the header of **every** read, not only those containing an N, so every N count is paired with the number of reads scanned on the same tile.
The report shows either:

- **N rate**: the fraction of reads with an N, for example at a given cycle or on a given tile. This is the default. Rates are comparable between tiles, lanes and samples with different cluster numbers.
- **N count**: the raw number of reads with an N.

Cycles are 1-based, so cycle 1 is the first base of the read.

### Sample, lane and read

- The **lane** comes from the read header.
- **Sample** and **read** (R1/R2/I1/...) come from the Illumina file name `<sample>_S<n>[_L<lane>]_<read>_001.fastq.gz`.
- Files that do not follow this convention use the file stem as the sample name and `R1` as the read.

## Interactive report

`nspector_report.html` embeds the Vega, Vega-Lite and vega-embed JavaScript, so it opens without a network connection (for example, copied from an HPC node).
The **Lane**, **Read** and **Metric** controls stay pinned at the top of the page.

| View | What it shows | Interaction |
| --- | --- | --- |
| Samples × lanes | N rate per sample and lane for the selected read. The selected lane is outlined. | Click a sample to select it. Shift-click adds more samples; double-click clears the selection. |
| N per cycle | One line per sample along the read. The selected samples are highlighted in red. | Drag to select a cycle window. |
| Flowcell tiles | Tiles laid out by surface, swath and segment, coloured by N bases per read **within the selected cycle window**. | Click a tile to select it. |
| N bases per read | Histogram (symlog) of N bases per read for the selected samples, or for all samples. | Follows the sample selection. |
| Tile × cycle | Heatmap of every tile at every cycle. Cells with no N are left empty. | Highlights the selected tile. |
| Spatial N density | N-containing reads binned over the X/Y coordinates of the selected tile (or of all tiles in the lane), with marginal histograms. | Follows the tile selection. |

The tile-level views (tile map, tile × cycle and spatial) sum over samples. Tile defects are physical, so they affect every sample on that lane.

Under the charts, the page lists:

- the tiles with the highest N rate;
- the cycles with the highest N rate;
- each input file with its read counts.

If the spatial data would exceed about 300,000 rows, only the tiles with the highest N rate keep a spatial map, and a warning is logged.

## Tables

Each table is written as `nspector_<table>.parquet`:

| Table | Grain | Columns |
| --- | --- | --- |
| `sample_lane` | sample, read, lane | `n`, `n_bases`, `total`, `n_rate` |
| `sample_cycle` | sample, read, lane, cycle (dense) | `n`, `total`, `n_rate` |
| `tile_summary` | sample, read, lane, tile | `surface`, `swath`, `segment`, `n`, `n_bases`, `total`, `n_rate`, `mean_n_per_affected_read` |
| `tile_cycle` | read, lane, tile, cycle (summed over samples; only cycles with an N) | `n`, `total`, `n_rate` |
| `n_per_read` | sample, read, lane, number of N bases | `n_reads` |
| `spatial` | read, lane, tile, `x_bin`, `y_bin` (summed over samples) | `n` |

In these tables:

- `n` is the number of reads with at least one N.
- `n_bases` is the number of N base calls.
- `total` is the number of reads scanned.

The raw per-file data (one row per N base, with tile, X, Y and cycle) is also kept as `<stem>_n_sequences.parquet`.

Example queries:

```python
import polars as pl

tiles = pl.read_parquet("nspector_results/nspector_tile_summary.parquet")
# worst tiles per lane, all samples combined
(tiles.group_by("lane", "tile")
      .agg(pl.col("n").sum(), pl.col("total").sum())
      .with_columns(rate=pl.col("n") / pl.col("total"))
      .sort("rate", descending=True)
      .head(10))
```

```sql
-- duckdb
SELECT read, lane, cycle, sum(n)::DOUBLE / sum(total) AS rate
FROM 'nspector_results/nspector_sample_cycle.parquet'
GROUP BY ALL ORDER BY rate DESC LIMIT 10;
```

## MultiQC

Nspector writes MultiQC custom content to `<output>/multiqc/`:

- `nspector_tile_rate_mqc.json`: a heatmap of % reads with N, with one row per sample, lane and read and one column per tile.
- `nspector_cycle_rate_mqc.json`: a line graph of % reads with N per cycle.
- `nspector_general_stats_mqc.tsv`: % reads with N and the N base count per sample, added to the General Statistics table.

```bash
multiqc nspector_results/multiqc
```

You can also point a larger MultiQC run at the output folder.

## Web interface

When nspector is run from `biomate web`, the result panel shows an **Open interactive report** button.
The server only serves reports produced by runs it started itself; each report is addressed by a random token, never by a file path.

## Output layout

```text
nspector_results/
├── nspector_report.html
├── nspector_sample_lane.parquet
├── nspector_sample_cycle.parquet
├── nspector_tile_summary.parquet
├── nspector_tile_cycle.parquet
├── nspector_n_per_read.parquet
├── nspector_spatial.parquet
├── SampleA_S1_L001_R1_001_n_sequences.parquet
├── ...
├── multiqc/
│   ├── nspector_tile_rate_mqc.json
│   ├── nspector_cycle_rate_mqc.json
│   └── nspector_general_stats_mqc.tsv
└── individual_tiles/                     # only with --png
    └── SampleA_S1_L001_R1_001_tile_1101.png
```

With `--png`, each input file also gets `<stem>_n_count_distribution.png` and `<stem>_n_distribution_per_cycle_and_tile.png`.

## Read header format

The module expects standard Illumina FASTQ read identifiers:

```text
@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filter>:<control>:<index>
```

Lane, tile, X and Y are taken from fields 3–6 (0-indexed) of the colon-separated identifier.

> _**Note:** A file is skipped with a warning when any read header deviates from this format. A file without any N bases is still included and contributes a rate of zero._
