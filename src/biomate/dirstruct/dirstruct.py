"""DirStruct main script."""

import argparse
import pathlib

from biomate.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
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
    parser_extract.add_argument(
        "--no-tags",
        action="store_true",
        help="Omit tags from the output file",
        default=False,
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
    """Main function."""
    if not args.output_file:
        args.quiet = True
    setup_logging(args)

    if args.command == "extract":
        # Delete the output file if it already exists
        if args.output_file and args.output_file.is_file():
            args.output_file.unlink()

        # Resolve the source path to ensure it is absolute
        args.source_path = args.source_path.resolve()
        # Walk through the directory structure
        for current_dir, _, files in pathlib.Path(args.source_path).walk():
            # Get the relative path of the directory to the source path
            current_dir = current_dir.relative_to(args.source_path.parent)
            # Case in which no output file is specified, just print to standard output
            if not args.output_file:
                if args.no_tags:
                    print(f"{current_dir}")
                    for file in files:
                        print(f"  {current_dir.joinpath(file)}")
                else:
                    print(f"{current_dir}\tdir")
                    for file in files:
                        print(f"  {current_dir.joinpath(file)}\tfile")
            # Write the directory structure to it
            else:
                # Create the parent directory if it does not exist
                if not args.output_file.parent.exists():
                    args.output_file.parent.mkdir(parents=True, exist_ok=True)
                # Write the directory structure to the output file
                with args.output_file.open("a") as output_file:
                    if args.no_tags:
                        output_file.write(f"{current_dir}\n")
                        for file in files:
                            output_file.write(f"  {current_dir.joinpath(file)}\n")
                    else:
                        output_file.write(f"{current_dir}\tdir\n")
                        for file in files:
                            output_file.write(f"{current_dir.joinpath(file)}\tfile\n")

    elif args.command == "create":
        if args.output_path:
            # Resolve the destination path to ensure it is absolute
            args.output_path = args.output_path.resolve()
            # Create the destination path if it does not exist
            if not args.output_path.exists():
                args.output_path.mkdir(parents=True, exist_ok=True)
        with args.source_file.open("r") as input_file:
            for line in input_file:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                # Parse line into path and optional type tag
                parts = line.split("\t")
                destination_path_str = parts[0].strip()
                destination_type_str = parts[-1].strip() if len(parts) > 1 else ""

                # Normalize destination type
                type_map = {"file": "file", "dir": "dir"}
                destination_type = type_map.get(destination_type_str, None)

                destination = args.output_path.joinpath(destination_path_str)

                # Create destination based on type
                try:
                    if destination_type == "file":
                        # Create parent directories and then the file
                        if not destination.parent.is_dir():
                            destination.parent.mkdir(parents=True, exist_ok=True)
                        destination.touch(exist_ok=True)
                    elif destination_type == "dir":
                        # Create the directory
                        destination.mkdir(parents=True, exist_ok=True)
                    else:
                        # No type specified: infer from file extension
                        # Files with an extension or special files (.DS_Store) are treated as files
                        if destination.suffix != "" or destination.name == ".DS_Store":
                            # Has extension or is .DS_Store - treat as file
                            if not destination.parent.is_dir():
                                destination.parent.mkdir(parents=True, exist_ok=True)
                            destination.touch(exist_ok=True)
                        else:
                            # No extension and not .DS_Store - treat as directory
                            destination.mkdir(parents=True, exist_ok=True)
                except OSError as e:
                    raise OSError(f"Failed to create '{destination}': {e}") from e
