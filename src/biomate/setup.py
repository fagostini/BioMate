import argparse
import logging


def setup_logging(args: argparse.Namespace) -> None:
    """
    Set up logging based on the verbosity level specified in the arguments.
    """
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    elif args.quiet:
        logging.basicConfig(level=logging.ERROR)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.debug(
        "Logging is set up with level: %s", logging.getLevelName(logging.root.level)
    )
