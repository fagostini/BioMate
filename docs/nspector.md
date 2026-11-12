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
    --output ./nspector_results
```

## Options

- `--input`: one or more FASTQ files (`.fastq`, `.fq`, optionally gzipped). Required.
- `--output`: output directory, created if missing (default: `./output`).
- `--engine`: `auto` (default), `python` or `rust`. `auto` uses the native `biomate_kernels` scanner when it is installed, otherwise Python. `rust` stops with an error when `biomate_kernels` is missing. Both engines produce identical results.
- `--threads`: Rust engine: files scanned in parallel (default `0` = one per available core). Python engine: extra decompression threads per file. The Rust engine also decompresses each gzip file on its own thread, ahead of parsing, and reads multi-member gzip (as written by BCL Convert).
- `--spatial-bin-size`: width of the spatial bins in read-name X/Y coordinate units (default: `2000`). Raise it if the report gets large.
- `--min-qual`: Phred score below which a called (non-N) base counts as low quality (default: `20`). See [Low-quality companion](#low-quality-companion).
- `--warn-rate`, `--fail-rate`: thresholds for the [anomaly flags](#anomaly-flags) (defaults: `0.01` and `0.05`).
- `--run-folder`: an Illumina run folder. Nspector reads per-tile cluster density, %PF and %≥Q30 from its `InterOp/` files and plots them against the N rate. Needs the optional extra: `pip install 'biomate[interop]'`.
- `--raw`: also write every N-containing read (lane, tile, X, Y, N cycles) as Parquet parts under `<stem>_n_sequences/`. This output grows with the number of reads and always uses the Python engine.
- `--png`: also write the static PNG charts for each input file.
- `--no-html`: skip the HTML report.
- `--no-multiqc`: skip the MultiQC files.

## Counts and rates

Nspector reads the header of **every** read, not only those containing an N, so every N count is paired with the number of reads scanned on the same tile.
The report shows either:

- **N rate**: the fraction of reads with an N, for example at a given cycle or on a given tile. This is the default. Rates are comparable between tiles, lanes and samples with different cluster numbers.
- **N count**: the raw number of reads with an N.
- **low-Q rate**: the fraction of reads with a low-quality called base (see below).

Cycles are 1-based, so cycle 1 is the first base of the read.

### Low-quality companion

In the same pass, nspector counts called bases (A, C, G or T, not N) whose Phred score is below `--min-qual`.
`lowq` columns hold the number of reads with such a base: at that cycle in the cycle tables, anywhere in the read in the tile and sample tables.
The two metrics separate two failure modes:

- **High N, normal low-Q**: the basecaller gave up on those clusters.
- **High low-Q, few N**: the signal degraded, but bases were still called.

### N-run shape

For each read with an N, nspector records the longest run of consecutive Ns, whether the read is all N, and where the first N falls: the first 10 cycles, the middle, or the last 10 cycles (the start wins on reads of 20 cycles or fewer).
Random Ns give runs of 1. Bubbles and dark cycles give long runs or all-N reads.

## Anomaly flags

Flags are computed from the tables after the scan, so they cost no extra reading. All thresholds use one scale: the fraction of reads with an N at a cycle.

| Kind | Raised when |
| --- | --- |
| `dark_cycle` | At least 80% of a lane's tiles are at or above `--warn-rate` at that cycle. |
| `tile` | The tile's mean per-cycle N rate, leaving out dark cycles, is at or above `--warn-rate`. |
| `tile_cycle` | One tile at one cycle is at or above `--warn-rate` and is a robust outlier among the lane's tiles at that cycle: `(rate - median) / (1.4826 × MAD) ≥ 5`. Cells already explained by a dark cycle, or by a tile flag at least as severe, are not listed again. |

A flag is `FAIL` when its rate is at or above `--fail-rate`, otherwise `WARN`.
The defaults come from a NextSeq 2000 run, where healthy cycles had about 0.1–0.15% of reads with an N and the problem cycles had 3–6%. Adjust them for your instruments.

Flags appear:

- in `nspector_flags.parquet`;
- at the top of the report, and as orange outlines on the tile map and the tile × cycle heatmap;
- in a status column of the worst-tiles table;
- in MultiQC, as an `N status` column (worst flag on the sample's lanes) and a flags table.

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
| N vs low-Q per cycle | Fraction of reads with an N, and with a low-quality called base, along the read. | Follows the lane, read and sample selection. |
| Longest N run, first N | Histogram of the longest N run per read, and where the first N falls. | Follows the lane, read and sample selection. |
| Tile × cycle | Heatmap of every tile at every cycle. Cells with no N and no low-Q base are left empty. Flagged cells and dark cycles are outlined in orange. | Highlights the selected tile. |
| N rate vs InterOp (with `--run-folder`) | One point per tile: N rate against %PF, cluster density or %≥Q30, coloured by lane. | Click a point to select that tile everywhere. |
| Spatial N density | N-containing reads binned over the X/Y coordinates of the selected tile (or of all tiles in the lane), with marginal histograms. | Follows the tile selection. |

The tile-level views (tile map, tile × cycle and spatial) sum over samples. Tile defects are physical, so they affect every sample on that lane.

Under the charts, the page lists:

- the tiles with the highest N rate;
- the cycles with the highest N rate;
- each input file with its read counts.

The embedded data depends on lanes, tiles and cycles, not on the number of reads. Two budgets keep it bounded:

- If the spatial data would exceed about 300,000 rows, only the tiles with the highest N rate keep a spatial map.
- If the tile × cycle data would exceed about 200,000 rows, cycles are merged into fixed windows in the report. The Parquet table keeps every cycle.

Both log a warning. With more than 200 samples, the cycle chart shows the median and interquartile range across samples, and draws only the samples you click. Charts render on a canvas, which stays fast with large heatmaps.

## Tables

Each table is written as `nspector_<table>.parquet`:

| Table | Grain | Columns |
| --- | --- | --- |
| `sample_lane` | sample, read, lane | `n`, `n_bases`, `lowq`, `total`, `n_rate`, `lowq_rate` |
| `sample_cycle` | sample, read, lane, cycle (dense) | `n`, `lowq`, `total`, `n_rate`, `lowq_rate` |
| `tile_summary` | sample, read, lane, tile | `surface`, `swath`, `segment`, `n`, `n_bases`, `lowq`, `total`, `n_rate`, `lowq_rate`, `mean_n_per_affected_read` |
| `tile_cycle` | read, lane, tile, cycle (summed over samples; only cycles with an N or a low-Q base) | `n`, `lowq`, `total`, `n_rate`, `lowq_rate` |
| `n_per_read` | sample, read, lane, number of N bases | `n_reads` |
| `spatial` | read, lane, tile, `x_bin`, `y_bin` (summed over samples; bin = coordinate // `--spatial-bin-size`) | `n` |
| `n_shape` | sample, read, lane, `longest_run`, `first_n`, `all_n` | `n_reads` |
| `flags` | one row per flag | `status`, `kind`, `read`, `lane`, `tile`, `cycle`, `rate`, `z`, `detail` |
| `tile_interop` (with `--run-folder`) | read, lane, tile (summed over samples) | `n`, `total`, `n_rate`, `density_k`, `pct_pf`, `pct_q30` |

In these tables:

- `n` is the number of reads with at least one N.
- `n_bases` is the number of N base calls.
- `lowq` is the number of reads with a called base below `--min-qual`.
- `total` is the number of reads scanned.

All tables are aggregated while the files are read, so memory use does not grow with the number of reads. Per-read data is only written with `--raw`:

```python
import polars as pl

raw = pl.read_parquet("nspector_results/SampleA_S1_L001_R1_001_n_sequences/*.parquet")
```

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
- `nspector_general_stats_mqc.tsv`: % reads with N, the N base count, % reads with a low-Q base and the worst flag status per sample, added to the General Statistics table.
- `nspector_flags_mqc.json`: a table of every flag (only written when there are flags).

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
├── nspector_n_shape.parquet
├── nspector_flags.parquet
├── nspector_tile_interop.parquet         # only with --run-folder
├── SampleA_S1_L001_R1_001_n_sequences/   # only with --raw
├── ...
├── multiqc/
│   ├── nspector_tile_rate_mqc.json
│   ├── nspector_cycle_rate_mqc.json
│   ├── nspector_general_stats_mqc.tsv
│   └── nspector_flags_mqc.json
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
