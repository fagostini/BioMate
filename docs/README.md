# Welcome to the BioMate

BioMate represents a collection of scripts and utility tools for everyday processing of biological data, conveniently gathered into a single Python package.

To check the package information and available sub-commands, you can run `biomate --help`.

## Available tools

- [Blabber](blabber.md) - FASTA, FASTQ or plain text sequence generator. It can be used to re-generate a demultiplexed flowcell folder using a sample sheet as input.
- [Dirstruct](dirstruct.md) - Extract or create a directory structure. It can be used to simulate nested folder structures for testing purposes.
- [Indexer](indexer.md) - Indexing tool for FASTA/FASTQ files. It can be used to search for expected indexes and get 'fuzzy' information.
- [FastRewind](fastrewind.md) - Convert a demultiplexed folder back into the Illumina NovaSeqXPlus output. It can be used to generate a subset for demultiplex testing.
- [Nspector](nspector.md) - Inspect FASTQ files for N-base content and their distribution across flowcell tiles and sequencing cycles.
- [Strainer](strainer.md) - Evaluate index mixing across flowcell lanes by identifying sample indexes from one lane appearing in undetermined reads of another.

## General options

In addition to the `--version` option, which prints the current version of the BioMate package, the following options can be provided before any sub-command (_e.g._ `biomate <OPTION> <SUBCOMMAND>`), and they will be propagated to the sub-command execution:

- `--verbose`: Increase verbosity of the output by setting the logging level to `DEBUG`. This will print additional information about the execution of the sub-command.
- `--quiet`: Decrease verbosity of the output by setting the logging level to `ERROR`. This will suppress most of the output, only printing errors.

> _**Note:** By default, the logging level is set to `INFO`, which will print general information about the execution of the sub-command, but not the detailed debug information._

> _**Important:** In some cases (e.g. when the output will be printed to the standard output), the `--quiet` option will be applied automatically to avoid cluttering the output with unnecessary information. In such cases, the logging level will be set to `ERROR` regardless of the `--verbose` or `--quiet` options._
