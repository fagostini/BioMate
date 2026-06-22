# Nspector

This module inspects one or more FASTQ files for N-base content and visualises how those bases are distributed across flowcell tiles and sequencing cycles.
It is useful for diagnosing quality issues tied to specific physical locations on the flowcell or to particular cycle positions.

## Usage

The minimal command is `biomate nspector --input <FASTQ_FILE>`.
Multiple input files can be provided at once; each is processed independently and produces its own set of output charts.

## Options

- `--input`: Path to one or more input FASTQ files to inspect. Gzipped (`.fastq.gz`) and uncompressed (`.fastq`) files are supported. This argument is required.
- `--output`: Path to the output directory where the charts will be saved. If the directory does not exist it will be created automatically (default: `./output`).

## Output

For each input file two chart images are written to the output directory, named after the input file stem:

- `<stem>_n_count_distribution.png` — a bar chart showing how many reads (y-axis, log scale) contain a given number of N bases (x-axis). Reads with no N bases are excluded.
- `<stem>_n_cycle_tile_distribution.png` — a line chart showing the number of N-containing reads per cycle (x-axis), with one coloured line per flowcell tile. This makes it easy to spot whether N bases cluster on a particular tile or emerge at a specific cycle.

> _**Note:** Files in which no reads contain an N base, or whose read headers do not conform to the standard Illumina format (`instrument:run:flowcell:lane:tile:x:y`), are skipped with a warning and do not produce output charts._

## Read Header Format

The module expects standard Illumina FASTQ read identifiers of the form:

```text
@<instrument>:<run>:<flowcell>:<lane>:<tile>:<x>:<y> <read>:<filter>:<control>:<index>
```

The tile, x and y coordinates are extracted from fields 5–7 (0-indexed) of the colon-separated identifier portion. Any FASTQ file whose headers deviate from this format will be skipped.

## Example

```bash
biomate nspector --input lane1_R1.fastq.gz lane2_R1.fastq.gz --output ./nspector_results
```

This will produce four files:

```text
nspector_results/
├── lane1_R1_n_count_distribution.png
├── lane1_R1_n_cycle_tile_distribution.png
├── lane2_R1_n_count_distribution.png
└── lane2_R1_n_cycle_tile_distribution.png
```
