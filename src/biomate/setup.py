import argparse
import logging


def setup_logging(args: argparse.Namespace) -> None:
    """
    Set up logging based on the verbosity level specified in the arguments.
    """
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    elif args.quiet:
        logging.getLogger().setLevel(logging.ERROR)
    else:
        logging.getLogger().setLevel(logging.INFO)

    logging.debug(
        "Logging is set up with level: %s", logging.getLevelName(logging.root.level)
    )
