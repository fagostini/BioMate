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

When no output is specified, the command will print the statistics for each matched pattern to the standard output. The output will look like this:

```text
0 CATGGGT 3902 0.39%%
0 CATGGTG 3584 0.36%%
0 CATGGTT 3471 0.35%%
0 CATGTTT 2717 0.27%%
0 CATGGGG 88 0.01%%
0 CATGTGT 82 0.01%%
0 CATGTGG 82 0.01%%
0 CATGTTG 79 0.01%%
```

The first column is the edit distance, the second column is the matched pattern, the third column is the number of matches, and the fourth column is the percentage of matches relative to the total number of sequences in the FASTQ file.
