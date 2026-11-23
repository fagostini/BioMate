# Triage

This module triages the undetermined reads of a demultiplexed flowcell.
It counts the header indexes of every undetermined FASTQ file **once**, then produces the reports of several modules from those counts, instead of reading the files once per module.

## Usage

The minimal command is `biomate triage --input-path <FLOWCELL_DIR>`.
The module expects to find the undetermined FASTQ files (`**/Undetermined*R1*.fastq.gz`) and a `SampleSheet.csv` inside the input directory.

## Options

The options are the same as for [Strainer](strainer.md): `--input-path`, `--output-path`, `--sample-sheet` and `--threads` (number of undetermined files counted in parallel).

## Output

All files are written to `--output-path`:

| File | Content |
|---|---|
| `<flowcell_name>_unexpected_indexes.csv`, `<flowcell_name>_lane_summary.csv` | The [Strainer](strainer.md) reports: sample indexes found in the undetermined reads of another lane, or reverse-complemented in any lane. Written under the same conditions as in Strainer. |
| `<flowcell_name>_barcode_findings.csv` | Undetermined barcodes explained by the sample sheet indexes of their lane (see [Barcode findings](#barcode-findings)). Always written, with only the header row when there is no finding. |
| `<flowcell_name>_top_undetermined.csv` | For each undetermined file, the 20 largest index groups, as in the [index discovery](indexer.md#index-discovery) with `--distance 1`. Columns: `File`, `Index` (the most frequent index of the group), `Count`, `Fraction` (of all the reads of the file) and `Variants` (number of distinct indexes in the group). Always written, with only the header row when there is no data. |
| `<flowcell_name>_always_undetermined.csv` | The [Waif](waif.md) report (default `Stats.json` mode). Only written when the input directory holds re-run datasets (`Stats.json`) sharing a lane; otherwise Waif is skipped with an info message. Datasets with different index lengths in a shared lane are compared on their common index length; a dataset with mixed index lengths within one lane stops the triage with an error, after the other reports are written. |

## Barcode findings

The 1000 most frequent indexes of each undetermined file (poly-G excluded) are compared with the indexes of the samples of the same lane. Indexes with a run of `GGGG` or `NNNN` in either half, and indexes equal to a sample's own barcode, are skipped. Each finding has a `Label`:

- `Index Mixing`: the i7 belongs to one sample and the i5 to another sample of the lane (`Sample_ID` and `Other_Sample_ID`).
- `Index Swap`: the i7 holds a sample's i5 and the i5 holds its i7, each as given or reverse-complemented.
- `Reverse` / `Complement`: one half is the reverse, or the complement, of a sample index. An index whose reverse or complement is itself is not compared.
- `Putative Sample`: no rule above matched, the lane holds a single project, and the count is within 1.5 standard deviations of the median sample count of the lane: likely a sample missing from the sheet.

Reverse-complemented indexes are not repeated here: the Strainer reports already cover them.

When the input directory holds exactly one `Reports/Demultiplex_Stats.csv` (BCL Convert), only indexes with more reads than the smallest sample of their lane are reported, and putative samples are searched. Without it, every match is reported and putative samples are skipped.

Columns: `Lane`, `Label`, `Undetermined_Index`, `Count`, `Sample_ID`, `Other_Sample_ID`, `Sample_Project` and `Matched` (which part matched what).

The rules come from flowceller's `check_barcodes`, which applies them to StatusDB documents; here they run on the FASTQ headers, before StatusDB has the flowcell.

## Notes

- The index groups are built from the 1000 most frequent indexes of each file, so the counts of a group leave out its rarest variants.
- Unlike Strainer, the index groups keep the poly-G indexes (`GGGG`): a large poly-G group is itself a useful signal.
- [Nspector](nspector.md) is not part of the triage: it reads the sample FASTQ files, not the undetermined ones.
