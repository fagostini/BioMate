"""BioMate main script."""

import argparse
import logging
import sys
from importlib.metadata import version

from rich.logging import RichHandler

from biomate import (
    blabber,
    dirstruct,
    fastrewind,
    index,
    nspector,
    strainer,
    triage,
    waif,
    web_interface,
)

try:
    __version__ = version(__name__)
except Exception:
    __version__ = "unknown"

__all__ = ["__version__", "main"]

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    handlers=[RichHandler(markup=True, rich_tracebacks=True)],
)


class CustomParser(argparse.ArgumentParser):
    """Parser class to enable help on error."""

    def error(self, message):
        """Show the error and help."""
        self.print_usage(sys.stderr)
        self.exit(2, f"error: {message}\n")


def _add_verbosity_flags(subparsers: argparse._SubParsersAction) -> None:
    """Accept --verbose/--quiet after any (nested) sub-command too."""
    for sub in subparsers.choices.values():
        # SUPPRESS: when omitted here, keep the value set before the sub-command
        sub.add_argument(
            "--verbose",
            action="store_true",
            default=argparse.SUPPRESS,
            help="Enable verbose output (debug level logging)",
        )
        sub.add_argument(
            "--quiet",
            action="store_true",
            default=argparse.SUPPRESS,
            help="Suppress output (error level logging only)",
        )
        for action in sub._actions:
            if isinstance(action, argparse._SubParsersAction):
                _add_verbosity_flags(action)


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
    nspector.nspector.init_parser(subparsers)
    strainer.strainer.init_parser(subparsers)
    triage.triage.init_parser(subparsers)
    waif.waif.init_parser(subparsers)
    web_interface.init_parser(subparsers)
    _add_verbosity_flags(subparsers)

    args = parser.parse_args()

    # Configure logging based on verbosity flags
    from biomate.setup import setup_logging

    if hasattr(args, "verbose") or hasattr(args, "quiet"):
        setup_logging(args)

    if not hasattr(args, "parse") or not hasattr(args, "run"):
        parser.print_help()
    else:
        try:
            args = args.parse(args)
        except argparse.ArgumentTypeError as e:
            # invalid input: show usage and the message, not a traceback
            parser.error(str(e))
        args.run(args)
