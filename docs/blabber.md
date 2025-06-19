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

```
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

When generating sequences based on a sample sheet, it is possible to simulate cases in which multiple recipes are present on the same lane. Running the demux process sequentially for each recipe will result in the unassigned files containing redundant reads belonging to the other recipes. By using the `taint` option, the unassigned reads contain ~10% of the reads from each sample present in the shared lane.
