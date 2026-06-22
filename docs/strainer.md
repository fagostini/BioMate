# Strainer

This module evaluates index mixing across flowcell lanes.
It scans the undetermined FASTQ files produced by a demultiplexing run and identifies reads whose index sequences match (within one mismatch) a sample that was assigned to a _different_ lane.
Such cross-lane contamination can occur when sample plates are accidentally loaded into the wrong lane position.

## Usage

The minimal command is `biomate strainer --input-path <FLOWCELL_DIR>`.
The module expects to find the undetermined FASTQ files and a `SampleSheet.csv` inside the input directory.

## Options

- `--input-path`: Path to the input flowcell directory. The directory must contain at least one file matching `**/Undetermined*R1*.fastq.gz`. This argument is required.
- `--output-path`: Path to the output directory where the results CSV will be written (default: current directory). The directory will be created automatically if it does not exist.
- `--sample-sheet`: Path to an explicit Sample Sheet file. If not provided, the module looks for `SampleSheet.csv` inside `--input-path`.
- `--threads`: Number of parallel worker processes used to read the undetermined FASTQ files. Default is `1` (single process). Increasing this value speeds up the extraction step on multi-core machines.

## Sample Sheet Format

The sample sheet must contain a `[Data]` or `[BCLConvert_Data]` section with at least the following columns:

| Column | Description |
|---|---|
| `Lane` | Lane number |
| `Sample_Project` | Project name |
| `index` | Index 1 sequence |
| `index2` | Index 2 sequence _(optional — single-index samples may omit this column)_ |

Lane values are normalised to the `LXXX` format (e.g. lane `1` becomes `L001`) so that they can be matched against the lane extracted from each undetermined filename.

## Matching Logic

For each sample in the sample sheet the module searches through the top-1000 most-frequent indexes from every _other_ lane's undetermined file.
An undetermined read is flagged as unexpected if:

- Its combined index string contains a sequence that matches `index` within **one substitution** (`s≤1`), **and**
- If `index2` is present, the combined index string also contains a sequence matching `index2` within one substitution.

Reads from the same lane as the sample are never flagged.

## Output

A single CSV file named `<flowcell_name>_unexpected_indexes.csv` is written to `--output-path`.
Each row describes one unexpected index observation:

| Column | Description |
|---|---|
| `Lane` | Lane in which the unexpected index was found |
| `Index` | The combined index string from the undetermined read |
| `Count` | Number of times this index was observed |
| `Sample_Project` | Project to which the matching sample belongs |
| `Lane_Project` | Lane to which the matching sample was assigned |
| `index1` | Index 1 of the matching sample |
| `index2` | Index 2 of the matching sample (empty for single-index samples) |

If no unexpected indexes are found the file is still written but will contain only the header row.

## Example

```bash
biomate strainer \
    --input-path /data/flowcells/20260310_LM43899_0385_A12GGASZR5 \
    --output-path ./strainer_results \
    --threads 4
```

The above command processes all undetermined FASTQ files in the flowcell directory using 4 parallel workers and writes the results to:

```text
strainer_results/
└── 20260310_LM43899_0385_A12GGASZR5_unexpected_indexes.csv
```

## Notes

- Only the top 1000 most-common undetermined indexes per lane are considered. Extremely rare cross-lane contamination events below this threshold will not be reported.
- Indexes containing four or more consecutive G bases (`GGGG`) are excluded from the analysis, as these are characteristic of sequencing errors on Illumina instruments rather than genuine index sequences.
- The module performs fuzzy matching with at most one substitution per index. Two substitutions or any insertion/deletion will not be matched.
