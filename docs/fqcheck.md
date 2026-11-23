# FqCheck

This module checks that FASTQ files are complete and that their mates hold the same reads.
Run it before a delivery, or after any step that rewrites FASTQ files (for example fastq-purge), to catch truncated or broken files that a checksum alone would not reveal.

## Usage

The minimal command is `biomate fqcheck --input-path <DIR>`.
The module searches `<DIR>` recursively for `*.fastq.gz` and `*.fastq` files. A link and its target count once, and a broken link is skipped with a warning.

## Options

- `--input-path`: directory to search (required).
- `--output-path`: directory where the report is written (default: current directory).
- `--md5`: also write `<dir_name>.md5`, an `md5sum`-compatible file with paths relative to `--input-path`. Check it with `cd <DIR> && md5sum -c <output>/<dir_name>.md5`.
- `--threads`: number of file groups checked in parallel (default: 1).

## Checks

Files are grouped by name with the read token (`_R1_`, `_R2_`, `_I1_`, `_I2_`) masked, so `A_S1_L001_R1_001.fastq.gz` and `A_S1_L001_R2_001.fastq.gz` form one group. A file without a read token is a group of its own.
Each group is read once, all mates in lockstep. A group has the status:

- `error` when a file is truncated or not valid gzip, a record is malformed (for example sequence and quality of different lengths), the mates hold different numbers of records, the read names of the mates are out of sync, or the group has no R1 file.
- `warning` when a file name does not follow the Illumina scheme (`_S<n>[_L<nnn>]_<R|I><1|2>_001.fastq[.gz]`) or the group holds no records.
- `ok` otherwise.

Only the first error of a group is reported.

## Output

`<dir_name>_fqcheck.csv`, one row per group, with the columns `Group`, `Files` (the mates, `;`-separated), `Records`, `Status` and `Message`.
The command exits with code 1 when any group has an error, after writing the report, so a pipeline can stop on it.
