"""BioMate main script."""

import argparse
import logging
import sys
from importlib.metadata import version

from rich.logging import RichHandler

from biomate import blabber, dirstruct, fastrewind, index

try:
    __version__ = version(__name__)
except Exception as e:
    raise e

__all__ = ["__version__", "biomate"]

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)


class CustomParser(argparse.ArgumentParser):
    """Parser class to enable help on error."""

    def error(self, message):
        """Show the error and help."""
        sys.stderr.write("\nError: %s\n\n\n" % message)
        self.print_help()
        print(
            "\nTo check the sub-commands help message use: biomate <sub-command> --help\n"
        )
        sys.exit(2)


def main():
    """Main function."""
    parser = CustomParser(
        description=f"BioMate {__version__}: A package for bioinformatics utilities.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose output (debug level logging)",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress output (error level logging only)",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
        help="Show the version of BioMate",
    )

    subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )
    blabber.blabber.init_parser(subparsers)
    dirstruct.dirstruct.init_parser(subparsers)
    index.index.init_parser(subparsers)
    fastrewind.fastrewind.init_parser(subparsers)

    args = parser.parse_args()
    if not hasattr(args, "parse") or not hasattr(args, "run"):
        parser.print_help()
    else:
        args = args.parse(args)
        args.run(args)
