# Blabber

This module allows to generate different flavours of FASTA, FASTQ or plain text sequences. It can be used to re-generate a flowcell folder using a sample sheet as input, or to simply generate random sequences for testing purposes.

## Usage

The minimal command is `biomate blabber`, which will generate 100 random 50 nucleotide long sequences in TEXT format (one sequence per line), using the standard nucleotide alphabet (A, C, G, T), and it will print them to the standard output.

## Options

- `--seq-number`: Number of sequences to generate. Default is 100.
- `--seq-length`: Length of each sequence. Default is 50.
- `--alphabet`: Alphabet to use for sequence generation. Default is 'ACGT'.
- `--output`: Path to the output file. If not specified, the output will be printed to the standard output. If a sample sheet is provided, this must be a directory where the flowcell folder will be created.
- `--format`: Format of the output sequences. Can be 'fasta', 'fastq', 'fastq-ext' or 'text'. Default is 'text'. The 'fastq-ext' format includes additional information in the FASTQ header.
- `--sample-sheet`: Path to a sample sheet file. If provided, the output will be generated based on the sample sheet, creating a flowcell folder structure. The sample sheet should be in the format used by Illumina sequencing systems, with the
- `--flowcell-id`: Full run folder name (`date_instrument_run_number_flowcell_id`) to use as the output folder; when given, the instrument, run number and flowcell in the read names are also pinned to its parts (default: a random identity is generated).
- `--taint`: Add reads from the samples to the Undetermined files of lanes with mixed sequencing designs (see below).
- `--taint-count`: Write exactly this many reads to the Undetermined files of each lane with an indexed design (implies `--taint`, deterministic selection, see below).
- `--random-seed`: Seed for the random number generator. Only affects the percentage-based `--taint` selection; `--taint-count` selection is seed-independent.

## Sample Sheet Usage

When a sample sheet is provided, the output will be generated based on the sample sheet, creating a flowcell folder structure. The sample sheet should be in the format used by Illumina sequencing systems, at least with the following columns:

- `FCID`: Flowcell ID.
- `Lane`: Lane number.
- `Sample_ID`: Sample ID.
- `Recipe`: Recipe name.
- `Sample_Project`: Sample project name.

> _**Note:** The `--seq-length` and `--format` options will be ignored when a sample sheet is provided, as the sequence lengths will be derived from the Recipe column in the sample sheet, and the format will be set to 'fastq-ext'._

The name of the output directory will have the following format: `<DATE>_<RAND_1>_<RAND_2>_A<FCID>`, where `<DATE>` is the current date in `YYYYMMDD` format, `<RAND_1>` is a random 2-character string followed by a random 5-digit number, and `<RAND_2>` is a random 4-digit number, and `<FCID>` is the Flowcell ID from the sample sheet.
Within the output flowcell folder, a `Demultiplexing` folder will be created, containing `Sample_Project` subfolders and the Undetermined FASTQ files for each lane. Each `Sample_Project` The `Sample_ID` subfolders that in turn contain the FASTQ files for each sample. The FASTQ files will be named according to the following pattern: `<Sample_Name>_<RAND_1>_S<RAND_2>_<Lane>_R[12]_001.fastq.gz`, where `<Sample_Name>` is derived from the corresponding column in the sample sheet, `<RAND_1>` is a 4-digit random number, `<RAND_2>` is the sample index (from 1 to the total number of samples in the flowcell), and `<Lane>` is the lane number from the sample sheet, converted to the format `LXXX` (e.g., `L001` for lane 1).

The `--seq-number` and `--alphabet` options will be used to generate the content of the FASTQ files, with the specified number of sequences and using the specified alphabet. The generated sequences will be random, but their identifiers will be consistent among paired reads (R1 and R2) to ensure they can be processed correctly in downstream applications.

The input sample sheet will be copied to the flowcell folder as `SampleSheet.csv`.

### Example Folder Structure

```text
20250613_ON46080_0992_A34YHADSA7
├── Demultiplexing
│   ├── A__Lpha_25_01
│   │   ├── Sample_P41583_101
│   │   │   ├── P41583_101_S1_L001_R1_001.fastq.gz
│   │   │   └── P41583_101_S1_L001_R2_001.fastq.gz
│   │   └── Sample_P41583_102
│   │       ├── P41583_102_S2_L001_R1_001.fastq.gz
│   │       └── P41583_102_S2_L001_R2_001.fastq.gz
│   └── Undetermined_S0_L001_R1_001.fastq.gz
│   └── Undetermined_S0_L001_R2_001.fastq.gz
└── SampleSheet.csv
```

### Simulating Mixed Design Recipes

When generating sequences based on a sample sheet, it is possible to simulate cases in which multiple recipes are present on the same lane. Running the demux process sequentially for each recipe will result in the unassigned files containing redundant reads belonging to the other recipes. Two options exist to populate the Undetermined files with such reads:

- `--taint` adds a percentage-based taint: roughly 10% of each sample's reads are flagged as candidates and, per lane, about 10% of the collected candidates (at least one) are written to the Undetermined files — around 1% of each sample's reads in total. The selection depends on `--random-seed`.
- `--taint-count N` (implies `--taint`) writes **exactly `N`** reads to the Undetermined files of every lane with an indexed design, spreading them round-robin over the lane's samples in sample-sheet order (bounded by the reads available in the lane, with a warning when fewer could be collected). Lanes whose samples all use a no-index design are skipped, because position-based demultiplexing cannot keep their reads undetermined. The selection is deterministic and independent of `--random-seed`, so the resulting counts are assertable in tests.

In both cases each tainted read is a full duplicate of a read already written to its sample's FASTQ files (same read name, sequence and qualities); with `--taint-count` the duplicate is placed at a unique tile/x/y position within the lane, so it never collides with the position of another read.

### Tile Distribution

When a sample sheet is provided, the available tiles (2 surfaces × 4 swaths × 16 tiles = 128) are distributed evenly across the lanes present in the sheet, so that each lane receives the same number of tiles (pool sizes differ by at most one when the tile count is not divisible by the number of lanes). Within each lane, the tile pool is shared among its samples — each sample's sequences are assigned to its tiles round-robin — so that every tile used in a lane appears in at least one sequence, regardless of how many samples or projects a lane contains.

All available tiles are used by default. If a lane's total number of sequences (`samples per lane × --seq-number`) is lower than its tile pool, the per-lane tile count is capped to the largest number of tiles that every lane can fully cover, so that each lane still uses exactly the same number of tiles and no used tile is left without sequences.

## Caveats

When simulating a flowcell run, there are some assumptions and tricks that do not really allow us to create a real-case dataset. These are mainly represented by the random numbers used in reads naming. Specifically, there is no check for duplicate entries (e.g. although unlikely, two reads might end up having the same tile, x and y positions). The balance of the tile distribution across lanes, however, is now guaranteed as described above.

When a tainted flowcell is re-converted with `fastrewind` and `bcl-convert`, undetermined reads belonging to lanes with shorter-than-full-length read designs may end with a run of `A` bases after the real read. This is an artifact of the BCL encoding: `fastrewind` writes the `N` padding cycles with the same base bits as `A` (the BCL base channel has no dedicated `N` slot), and `bcl-convert` trims undetermined clusters differently from sample-sheet reads, which are trimmed to their planned read length. Sample files are unaffected, and the undetermined read count per lane is exact.
