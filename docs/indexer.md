# Index

This module allows you to scan a fastq file, either gzipped or uncompressed, and extract the sequences that match a given pattern, or a list of patterns.

## Usage

The command `biomate index` takes at least two arguments: the path to a FASTQ file (`--input`) and the pattern to search for. The pattern can be a single string or a regular expression (`--index-regex`), or a path to a file containing a list of patterns or regular expressions (`--index-file`). If no output (`--output`) is specified, only the statistics for each matched pattern will be printed to the standard output; otherwise, both the statistics for each matched pattern (`pattern<INDEX>_matches.txt`) and for each edit distance (`pattern<INDEX>_errors.txt`) will be written to files, where `<INDEX>` is the index of the pattern in the input file or `1` if a single pattern is provided.

## Options

- `--input`: Path to the input FASTQ file. It can be gzipped or uncompressed.
- `--output`: Path to the output directory where the results will be saved. If not specified, the results will be printed to the standard output. The output directory will contain two files: `indexing_patterns.txt` with the statistics for each matched pattern, and `indexing_errors.txt` with the statistics for each edit distance.
- `--index-regex`: Regular expression to search for in the sequences.
- `--index-file`: Path to a file containing a list of patterns or regular expressions to search for in the sequences. Each line in the file should contain a single pattern or regular expression.
- `--distance`: Maximum edit distance to consider when matching patterns. Default is 0 (exact match).
- `--error-type`: Type of error to consider when calculating the edit distance. Can be insertions (`i`), deletions (`d`), substitutions (`s`), or any error (`e`). Default is `s`.

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
