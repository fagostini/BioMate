"""Logging setup."""

import argparse
import logging


def setup_logging(args: argparse.Namespace) -> None:
    """
    Set up logging based on verbosity arguments.

    Configures the root logger's level based on the --verbose and --quiet flags.
    Verbose takes precedence over quiet. Ensures at least one handler exists for
    log output (falls back to StreamHandler if none configured).

    Args:
        args: Parsed command-line arguments with 'verbose' and 'quiet' attributes

    Raises:
        AttributeError: If args is missing 'verbose' or 'quiet' attributes
    """
    # Validate required attributes exist
    if not hasattr(args, "verbose") or not hasattr(args, "quiet"):
        raise AttributeError(
            "Arguments must have 'verbose' and 'quiet' boolean attributes"
        )

    # Determine log level (verbose takes precedence over quiet)
    if args.verbose:
        log_level = logging.DEBUG
    elif args.quiet:
        log_level = logging.ERROR
    else:
        log_level = logging.INFO

    root_logger = logging.getLogger()
    root_logger.setLevel(log_level)

    # Ensure at least one handler exists for output
    if not root_logger.handlers:
        handler = logging.StreamHandler()
        handler.setLevel(log_level)
        root_logger.addHandler(handler)

    logging.debug(
        "Logging configured with level: %s",
        logging.getLevelName(root_logger.level),
    )
