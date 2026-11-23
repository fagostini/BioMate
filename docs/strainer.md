# Strainer

> _**Deprecated:** `biomate strainer` now runs [Triage](triage.md), which writes the same reports plus its own, and will be removed in the next minor release. The search described below is unchanged._

This module evaluates index mixing across flowcell lanes.
It scans the undetermined FASTQ files produced by a demultiplexing run and identifies reads whose index sequences match (within one mismatch) a sample that was assigned to a _different_ lane.
Such cross-lane contamination can occur when sample plates are accidentally loaded into the wrong lane position.

## Usage

The minimal command is `biomate strainer --input-path <FLOWCELL_DIR>`.
The module expects to find the undetermined FASTQ files and a `SampleSheet.csv` inside the input directory.

## Options

- `--input-path`: Path to the input flowcell directory. The directory must contain at least one file matching `**/Undetermined*R1*.fastq.gz`, and every such file must be split by lane, with the lane in its name (`*_LXXX_*`, e.g. `Undetermined_S0_L001_R1_001.fastq.gz`); files without a lane (BCL Convert `--no-lane-splitting`) are rejected. This argument is required.
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

For each sample in the sample sheet the module searches through the 1000 most-frequent indexes of each undetermined file of every _other_ lane.
An undetermined read is flagged as unexpected if:

- Its combined index string contains a sequence that matches `index` within **one substitution** (`s≤1`), **and**
- If `index2` is present, the combined index string also contains a sequence matching `index2` within one substitution.

With the indexes as written in the sample sheet (`Orientation` = `expected`), reads from the same lane as the sample are never flagged. A hit is also dropped when the same index pair is assigned to a sample in the lane where it was found, since it is then expected there.

### Index orientation

Each sample is also searched with its indexes reverse-complemented: `index2_rc`, `index1_rc` and `both_rc` (only `index1_rc` for single-index samples).
These are searched in every lane, the sample's own included: an index entered in the wrong orientation leaves the sample undetermined in its own lane.
A reverse-complemented index must match exactly (no substitution) and in its own half of the combined string (`index` before the `+`, `index2` after it): a wrongly oriented sample shows up as its exact reverse complement, in bulk.
A swap of `index` and `index2` is not reported separately: with the indexes as written in the sample sheet, the matching ignores the position of each index within the combined string.

## Output

Two CSV files are written to `--output-path`. The first, `<flowcell_name>_unexpected_indexes.csv`, lists the matches.
Each row describes one unexpected index observation:

| Column | Description |
|---|---|
| `Lane` | Lane in which the unexpected index was found |
| `Undetermined_Index` | The combined index string from the undetermined read |
| `Count` | Number of times this index was observed |
| `Fraction` | `Count` divided by the undetermined reads searched (the reads of the top 1000 non-poly-G indexes of each file): those of the other lanes for `expected`, those of all lanes for the reverse-complement orientations |
| `Sample_Project` | Project to which the matching sample belongs |
| `Lane_Project` | Lane to which the matching sample was assigned |
| `index1` | Index 1 of the matching sample |
| `index2` | Index 2 of the matching sample (empty for single-index samples) |
| `Orientation` | `expected`, or the reverse-complement combination that matched (`index2_rc`, `index1_rc`, `both_rc`) |

The second, `<flowcell_name>_lane_summary.csv`, has one row per lane and orientation, with the columns `Lane`, `Lane_Count`, `Lane_Fraction` (sums over the distinct undetermined indexes) and `Orientation`.

If the search finds no candidate at all, neither file is written. If it finds candidates but all of them are then dropped as expected, both files are written with only their header row.

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
├── 20260310_LM43899_0385_A12GGASZR5_lane_summary.csv
└── 20260310_LM43899_0385_A12GGASZR5_unexpected_indexes.csv
```

## Notes

- To run Strainer together with the index discovery and Waif from a single pass over the undetermined files, see [Triage](triage.md).
- Only the 1000 most-common undetermined indexes of each file are considered. Extremely rare cross-lane contamination events below this threshold will not be reported.
- Indexes containing four or more consecutive G bases (`GGGG`) are excluded from the analysis, as these are characteristic of sequencing errors on Illumina instruments rather than genuine index sequences.
- The module performs fuzzy matching with at most one substitution per index. Two substitutions or any insertion/deletion will not be matched.
