# Waif

This module reports the barcodes that remain **undetermined in every demultiplexing run of a lane**.
When a flowcell lane is demultiplexed multiple times (for example, successive bcl-convert re-runs with an updated sample sheet), barcodes that are undetermined in *all* of those runs are strong candidates for index hopping or index mis-assignment.
For each such barcode the module reports an estimated read count: the mean of the per-dataset counts, rounded up to a multiple of 10.

## Usage

The minimal command is `biomate waif --input-path <DIR>`.
By default the module reads the `UnknownBarcodes` field of every `Stats.json` file found under the input directory.

## Options

- `--input-path`: Path to the input directory. In default mode it must contain at least one `<Dataset>/Stats/Stats.json` file; in exact mode it must contain the undetermined FASTQ files. This argument is required.
- `--output-file`: File to write the results CSV to. If not specified, the results are written to stdout (and logging is suppressed automatically). The parent directory is created automatically if it does not exist.
- `--exact`: Enable exact counting of undetermined indexes by reading the undetermined FASTQ files instead of relying on `Stats.json`. Requires a Sample Sheet to identify the lanes that carry samples from two or more projects.
- `--sample-sheet`: Path to an explicit Sample Sheet file (only used in exact mode). If not provided, the module looks for `SampleSheet.csv` inside `--input-path`.
- `--block-cache DIR`: Directory for a small gzip member-offset index of the undetermined files; enables forked per-file workers that decompress and count in parallel (effective when there are fewer files than cores). Each file is scanned once to locate member boundaries; the index file is tiny and the source is never copied or re-compressed. Only supported on Linux.

## Input Format

### Default mode (Stats.json)

The module globs `**/Stats.json` under `--input-path` and reads the `UnknownBarcodes` field, which is expected to look like:

```json
{
    "UnknownBarcodes": [
        {"Lane": 1, "Barcodes": {"AAAA+AAAA": 200, "CCCC+CCCC": 100}}
    ]
}
```

The dataset name is the directory that contains the `Stats` folder (i.e. `input_path/<Dataset>/Stats/Stats.json`).
Files below any `Demultiplexing/` directory are excluded, as that is the primary bcl-convert demux output rather than a re-run.
`Stats.json` files without an `UnknownBarcodes` field are skipped with a warning.

### Exact mode

The module globs `**/Undetermined_*_L00<Lane>_R1_*.fastq.gz` under `--input-path` and counts the index pairs in the read headers (only R1 files are used, as the index information is in the R1 headers).
Files below any `Demultiplexing/` directory are excluded.
If bcl-convert split a dataset into several 1M-read parts (`_001`, `_002`, ...), their counts are summed, so each (Lane, Dataset, Barcode) combination is unique.

The Sample Sheet must be in either the legacy `FCID,...` format or the modern `Lane,...` format (the one used by `[BCLConvert_Data]` sheets).
Rows without any index and rows with an empty `Sample_Project` are ignored; a lane is selected when its rows carry samples from two or more distinct projects.

## Output

A CSV with the following columns is written to `--output-file` (or to stdout):

| Column | Description |
|---|---|
| `Lane` | Lane number |
| `Barcode` | The undetermined barcode (`index1` or `index1+index2`) |
| `Count` | Estimated read count: mean of the per-dataset counts, rounded up to a multiple of 10 |
| `Total %` | Percentage of the lane's total reported count |

Rows are sorted by lane, then by count (descending), then by barcode.
If no barcode is undetermined in every run of a lane (or no lane is shared), the file is still written but contains only the header row.

## Counting Strategies (Exact Mode)

Each undetermined FASTQ file is counted with the fastest strategy that is correct for it:

- **Serial**: one fused streaming pass (decompress + count) for small files.
- **Threaded**: the `zlib` inflate C call (which releases the GIL) runs in a producer thread overlapped with the counting pass, for large files.
- **Forked (with `--block-cache`)**: a small index of the gzip member byte offsets is built and cached (validated against the source by path, mtime and size); the file is then counted by one forked worker per contiguous run of members. This is only used when every member boundary falls on a FASTQ line boundary; otherwise the module falls back to the serial pass. A multi-member file whose boundary falls mid-header is never split.

Per-file counting runs in a process pool (one worker per file, up to the core count) and the output order matches the file order, so the results are deterministic.

## Example

```bash
biomate waif \
    --input-path /data/flowcells/20260310_LM43899_0385_A12GGASZR5 \
    --output-file ./waif_results/undetermined.csv

biomate waif \
    --input-path /data/flowcells/20260310_LM43899_0385_A12GGASZR5 \
    --exact \
    --sample-sheet /data/sample_sheets/SampleSheet.csv \
    --block-cache /tmp/waif_cache
```

The first command processes all `Stats.json` files under the flowcell directory and writes the results to `undetermined.csv`.
The second command counts the undetermined FASTQ headers exactly, using the given sample sheet to select the lanes, and caches the gzip member-offset index in `/tmp/waif_cache`.

## Notes

- The reported count is an estimate: it is the mean of the per-dataset counts rounded up to a multiple of 10, with a single rounding applied at the end.
- Datasets may be dual-index (`index1+index2`) or single-index (`index1` only). When every dataset of a lane is dual-index the intersection is computed on both indexes; as soon as a single-index dataset joins, the comparison drops to index1 only (the longer index1 is trimmed to the shorter) and the result is reported in the single-index form.
- A (Dataset, Lane) must be uniformly single- or dual-index; mixing the two within one dataset is a data error and aborts the run.
- Two datasets with equal total index length but crossed per-index lengths (e.g. 10+8 vs 8+10) can leave one index untruncated, in which case that lane's intersection comes out empty.
- `--block-cache` is only supported on Linux (it relies on `os.fork`); on other platforms the argument is rejected.
- The web interface runs each tool with a 10-minute timeout; large exact-mode runs should be done from the command line.