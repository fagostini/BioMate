# Welcome to the BioMate

BioMate represents a collection of scripts and utility tools for everyday processing of biological data, conveniently gathered into a single Python package.

To check the package information and available sub-commands, you can run `biomate --help`.

## Available tools

- [Blabber](blabber.md) - FASTA, FASTQ or plain text sequence generator. It can be used to re-generate a flowcell folder using a sample sheet as input.
- [Dirstruct](dirstruct.md) - Extract or create a directory structure. It can be used to simulate nested folder structures for testing purposes.
- [Index](index.md) - Indexing tool for FASTA/FASTQ files. It can be used to search for expected indexes and get 'fuzzy' information.

## General options

In addition to the `--version` option, which prints the current version of the BioMate package, the following options can be provided before any sub-command (_e.g._ `biomate <OPTION> <SUBCOMMAND>`), and they will be propagated to the sub-command execution:

- `--verbose`: Increase verbosity of the output by setting the logging level to `DEBUG`. This will print additional information about the execution of the sub-command.
- `--quiet`: Decrease verbosity of the output by setting the logging level to `WARNING`. This will suppress most of the output, only printing warnings and errors.

> _**Note:** By default, the logging level is set to `INFO`, which will print general information about the execution of the sub-command, but not the detailed debug information._

> _**Important:** In some cases (e.g. when the output will be printed to the standard output), the `--quiet` option will be applied automatically to avoid cluttering the output with unnecessary information. In such cases, the logging level will be set to `WARNING` regardless of the `--verbose` or `--quiet` options._
