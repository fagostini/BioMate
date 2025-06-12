import argparse
import logging
import sys
from importlib.metadata import version

from rich.logging import RichHandler

from biomate import blabber, dirstruct, index

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
    def error(self, message):
        sys.stderr.write("error: %s\n" % message)
        self.print_help()
        sys.exit(2)


def main():
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

    subparsers = parser.add_subparsers(
        title="sub-commands",
        help="Access the help page for a sub-command with: sub-command -h",
    )
    blabber.blabber.init_parser(subparsers)
    dirstruct.dirstruct.init_parser(subparsers)
    index.indexing.init_parser(subparsers)

    args = parser.parse_args()
    if not hasattr(args, "parse") or not hasattr(args, "run"):
        parser.print_help()
    else:
        args = args.parse(args)
        args.run(args)
