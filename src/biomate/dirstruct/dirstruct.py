import argparse
import pathlib

from biomate.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Indexing tool for FASTA/FASTQ files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Indexing tool for FASTA/FASTQ files",
    )
    subcommands = parser.add_subparsers(dest="command", required=True)
    parser_extract = subcommands.add_parser(
        "extract", help="Extract the directory structure from a path"
    )
    parser_extract.add_argument(
        "--source-path",
        type=pathlib.Path,
        help="Path to the directory to index",
        required=True,
    )
    parser_extract.add_argument(
        "--output-file",
        type=pathlib.Path,
        help="Path to the output file where the directory structure will be saved",
        required=False,
    )
    parser_create = subcommands.add_parser(
        "create", help="Create the directory structure from a file"
    )
    parser_create.add_argument(
        "--source-file",
        type=pathlib.Path,
        help="Path to the file containing the directory structure",
        required=True,
    )
    parser_create.add_argument(
        "--output-path",
        type=pathlib.Path,
        help="Path to the destination directory where the structure will be created",
        default=pathlib.Path("./"),
        required=False,
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if args.command == "extract":
        if not args.source_path.is_dir():
            raise argparse.ArgumentTypeError(
                f"The source path '{args.source_path}' does not exist or is not a directory."
            )
    else:
        if not args.source_file.is_file():
            raise argparse.ArgumentTypeError(
                f"The source file '{args.source_file}' does not exist or is not a file."
            )
    return args


def main(args: argparse.Namespace) -> None:
    setup_logging(args)
    if args.command == "extract":
        # Delete the output file if it already exists
        if args.output_file and args.output_file.is_file():
            args.output_file.unlink()

        # Resolve the source path to ensure it is absolute
        args.source_path = args.source_path.resolve()
        # Walk through the directory structure
        for dir, _, files in pathlib.Path(args.source_path).walk():
            # Get the relative path of the directory to the source path
            dir = dir.relative_to(args.source_path.parent)
            # Case in which no output file is specified, just print to standard output
            if not args.output_file:
                print(f"{dir}")
                for file in files:
                    print(f"  {dir.joinpath(file)}")
            # Write the directory structure to it
            else:
                # Create the parent directory if it does not exist
                if not args.output_file.parent.exists():
                    args.output_file.parent.mkdir(parents=True, exist_ok=True)
                # Write the directory structure to the output file
                with args.output_file.open("a") as output_file:
                    output_file.write(f"{dir}\n")
                    for file in files:
                        output_file.write(f"{dir.joinpath(file)}\n")

    elif args.command == "create":
        if args.output_path:
            # Resolve the destination path to ensure it is absolute
            args.output_path = args.output_path.resolve()
            # Create the destination path if it does not exist
            if not args.output_path.exists():
                args.output_path.mkdir(parents=True, exist_ok=True)
        with args.source_file.open("r") as input_file:
            for line in input_file:
                destination = args.output_path.joinpath(line.strip())
                # Warning: This will not correctly handle files with no extension
                if destination.suffix == "":
                    destination.mkdir(parents=True, exist_ok=True)
                else:
                    destination.touch(exist_ok=True)
