# Index

This module allows you to scan one or more FASTQ files, either gzipped or uncompressed, and extract the sequences that match a given pattern, or a list of patterns. Without a pattern, it counts and groups all the indexes found in the read headers.

## Usage

The command `biomate index` requires the path to one or more FASTQ files (`--input`). With a pattern it searches for it; without one it runs [index discovery](#index-discovery). The pattern can be a single string or a regular expression (`--index-regex`), or a path to a file containing a list of patterns or regular expressions (`--index-file`). If no output (`--output`) is specified, only the statistics for each matched pattern will be printed to the standard output; otherwise, both the statistics for each matched pattern (`pattern<INDEX>_matches.txt`) and for each edit distance (`pattern<INDEX>_errors.txt`) will be written to files, where `<INDEX>` is the index of the pattern in the input file or `1` if a single pattern is provided.

## Options

- `--input`: Path(s) to the input FASTQ file(s). They can be gzipped or uncompressed. Results are kept per file: with several inputs, standard output prints a `# <path>` line before each file's results, and `--output` gets one subdirectory per input file, named after the file without its `.fastq`/`.fq`/`.fasta`/`.fa` suffix (optionally followed by `.gz`); input file names must therefore be distinct.
- `--output`: Path to the output directory where the results will be saved. If not specified, the results will be printed to the standard output. With a pattern, the output directory contains `pattern<INDEX>_matches.txt` (statistics for each matched sequence) and `pattern<INDEX>_errors.txt` (statistics for each edit distance) per pattern; with index discovery, it contains `indexes_matches.txt` and `indexes_groups.txt`.
- `--index-regex`: Pattern or regular expression to search for in the sequences. Patterns are uppercased, and a `+` separates the two parts of a dual index. Cannot be combined with `--index-file`.
- `--index-file`: Path to a file containing a list of patterns or regular expressions to search for in the sequences. Each line in the file should contain a single pattern or regular expression.
- `--engine`: Matching engine, `regex` (default) or `rust`. `rust` needs the optional `biomate_kernels` package and falls back to `regex` when the package is missing or the pattern uses regular-expression syntax. Ignored by index discovery.
- `--distance`: Maximum edit distance to consider when matching patterns. Default is 0 (exact match).
- `--min-count`: Index discovery only. Indexes seen fewer times are not grouped or listed; their read and index totals are reported in one `# ... below --min-count ...` line (at the end of standard output, or of `indexes_groups.txt`). Default is 1 (keep all). On large files most distinct indexes are sequencing noise seen once or twice, so a small cutoff (e.g. 10) cuts grouping time, memory and output size by orders of magnitude.
- `--threads`: Number of input files processed in parallel. Only used with several inputs and `--output`; standard output is always written one file at a time. Default is 1.
- `--error-type`: Type of error to consider when calculating the edit distance. Can be insertions (`i`), deletions (`d`), substitutions (`s`), or any error (`e`). Default is `s`.

## Index discovery

When neither `--index-regex` nor `--index-file` is given, the module reads the index from each read header (the text after the last `:` of the comment, e.g. `ACGTACGT+TTTTGGGG` in `@... 1:N:0:ACGTACGT+TTTTGGGG`) and groups similar indexes. Indexes are visited from most to least frequent: each joins the first, more frequent, index of the same length within `--distance` substitutions (Hamming distance; `--error-type` and `--engine` are ignored), otherwise it starts a new group. With `--distance 0` this is a plain count of identical indexes.

Index discovery reads FASTQ input only (gzipped or uncompressed, four lines per record); pattern search also accepts FASTA. The header format is checked on the first record: if it carries no index, a warning is logged and no index is reported for that file. Percentages are relative to the number of reads in the file.

Each output line is `<group seed> <distance> <index> <count> <percentage>`, where `<group seed>` is the sequence of the most frequent index of the group; groups are sorted by their total count. With `--output`, these lines go to `indexes_matches.txt`, and `indexes_groups.txt` holds one `<group seed> <total count> <percentage>` line per group.

## Example of Standard Output

When no output is specified, the command will print the statistics for each matched pattern to the standard output. As an example, the output of the command `biomate index --input Test_1M.fastq.gz --index-regex "AACGTGAT[CG]TG" --distance 1` would look like this:

```text
0 AACGTGATGTG 40 0.00%
1 AAGGTGATGTG 3020 0.30%
1 CACGTGATGTG 105 0.01%
1 ATCGTGATGTG 89 0.01%
1 AACGTGAAGTG 78 0.01%
1 AACGTGATGTT 60 0.01%
1 ACCGTGATGTG 37 0.00%
1 AATGTGATGTG 36 0.00%
1 AACGTGAACTG 24 0.00%
1 AGCGTGATGTG 24 0.00%
1 AAAGTGATGTG 16 0.00%
1 AACGTGGTGTG 7 0.00%
1 AACGTGGTCTG 5 0.00%
1 GACGTGATGTG 5 0.00%
1 AACGCGATGTG 5 0.00%
1 AACGAGATGTG 4 0.00%
1 AACGTGAGGTG 2 0.00%
1 AACTTGATGTG 1 0.00%
1 AACGTTATCTG 1 0.00%
1 AACGTGATCTC 1 0.00%
1 AACGCGATCTG 1 0.00%
1 TACGTGATGTG 1 0.00%
1 AACTTGATCTG 1 0.00%
```

The first column is the edit distance, the second column is the matched pattern, the third column is the number of matches, and the fourth column is the percentage of matches relative to the total number of sequences in the FASTQ file.
