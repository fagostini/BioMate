"""DirStruct main script."""

import argparse
import logging
import os
import pathlib

from biomate.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Extract or re-create a directory structure",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Extract or re-create a directory structure",
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
    # extract without --output-file prints the tree to stdout: keep logs out of it
    if args.command == "extract" and not args.output_file:
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
        n_created = n_existing = n_conflict = 0
        with args.source_file.open("r") as input_file:
            for line in input_file:
                line = line.strip()
                if not line:
                    continue  # Skip empty lines

                # Parse line into path and optional type tag
                parts = line.split("\t")
                destination_path_str = parts[0].strip()
                destination_type = parts[-1].strip() if len(parts) > 1 else ""

                destination = args.output_path.joinpath(destination_path_str)
                # normpath, not resolve: existing symlinks inside the output stay "existing"
                if not pathlib.Path(os.path.normpath(destination)).is_relative_to(
                    args.output_path
                ):
                    n_conflict += 1
                    logging.warning(
                        f"'{destination_path_str}' is outside the output path; skipping"
                    )
                    continue

                if destination_type in ("file", "dir"):
                    is_file = destination_type == "file"
                else:
                    # No type specified: infer from file extension
                    # Files with an extension or special files (.DS_Store) are treated as files
                    is_file = (
                        destination.suffix != "" or destination.name == ".DS_Store"
                    )
                kind = "file" if is_file else "directory"

                # Leave pre-existing paths untouched (symlinks included, even broken ones)
                if destination.exists() or destination.is_symlink():
                    if destination.is_dir() != is_file:
                        n_existing += 1
                        logging.debug(f"Keeping existing {kind} '{destination}'")
                    else:
                        n_conflict += 1
                        logging.warning(
                            f"'{destination}' exists but is not a {kind}; skipping"
                        )
                    continue

                try:
                    if is_file:
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        destination.touch()
                    else:
                        destination.mkdir(parents=True)
                except (FileExistsError, NotADirectoryError) as e:
                    # A parent component exists as a file
                    n_conflict += 1
                    logging.warning(
                        f"Cannot create {kind} '{destination}': {e}; skipping"
                    )
                    continue
                except OSError as e:
                    raise OSError(f"Failed to create '{destination}': {e}") from e
                n_created += 1

        logging.info(
            f"Created {n_created} entries, kept {n_existing} existing, "
            f"skipped {n_conflict} conflicting"
        )
