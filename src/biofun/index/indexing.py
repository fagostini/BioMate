import argparse
import logging
import pathlib
from collections import Counter, defaultdict

import dnaio
import regex as re

from biofun.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Indexing tool for FASTA/FASTQ files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Indexing tool for FASTA/FASTQ files",
    )
    parser.add_argument(
        "--input",
        type=pathlib.Path,
        required=True,
        help="Path to the input FastQ file",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        required=False,
        default=None,
        help="Path to the output path where the results will be saved",
    )
    index_group = parser.add_mutually_exclusive_group(required=True)
    index_group.add_argument(
        "--index-regex",
        type=str,
        default=None,
        help="Regular expression to search in the FastQ file",
    )
    index_group.add_argument(
        "--index-file",
        type=pathlib.Path,
        default=None,
        help="Path to the index file. One index per line",
    )
    distance_group = parser.add_argument_group("Distance Options")
    distance_group.add_argument(
        "--distance",
        type=int,
        default=1,
        help="Maximum distance between subject and query for sequences",
    )
    distance_group.add_argument(
        "--distance-type",
        type=str,
        choices=["i", "d", "s", "e"],
        default="s",
        help="Type of distance to use: 'i' for insertion, 'd' for deletion, 's' for substitution, 'e' for any edit",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s 0.1.0",
        help="Show the version of the indexing tool",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if not args.input.is_file():
        raise FileNotFoundError(f"Input file {args.input} does not exist.")
    return args


def compile_index_regex(index_regex, distance, error_type):
    """
    Compile the index regex with the specified distance and error type.
    """
    if error_type == "e":
        error_pattern = f"e<={distance}"
    else:
        error_pattern = ",".join(
            [
                f"{x}<={distance}" if x in error_type else f"{x}<=0"
                for x in ["i", "d", "s"]
            ]
        )

    split_regex = index_regex.upper().split("+")
    index_regex = ".*".join([f"({x}){{{error_pattern}}}" for x in split_regex])
    return re.compile(f"({index_regex}){{{error_pattern}}}")


def write_results(results, output_path):
    """Write the results to the specified output path."""
    if not output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)
    max_len = max(
        [len(" ".join(list(seq))) for _, value in results.items() for seq in value]
    )
    with open(output_path.joinpath("indexing_sequences.txt"), "w") as output_file:
        for key, value in results.items():
            for seq, cnt in value.most_common():
                output_file.write(f"{key} {' '.join(list(seq)):<{max_len}} {cnt}\n")

    with open(output_path.joinpath("indexing_errors.txt"), "w") as output_file:
        for key, value in results.items():
            counts = sum(value.values())
            output_file.write(f"{key} {counts}\n")


def main(args: argparse.Namespace) -> None:
    if not args.output:
        args.quiet = False
        args.quiet = True
    setup_logging(args)

    logging.info(f"Input file: {args.input}")
    input_type = (
        "FASTQ"
        if re.match(r".*\.f(ast)?q(.gz)$", args.input.name)
        else "FASTA"
        if re.match(r".*\.f(ast)?a(.gz)?$", args.input.name)
        else None
    )
    logging.info(f"Output type: {input_type}")
    logging.info(f"Output file: {args.output}")
    if args.index_regex:
        logging.info(f"Using regex for indexing: {args.index_regex}")
    if args.index_file:
        logging.info(f"Using index file: {args.index_file}")

    regex_pattern = compile_index_regex(
        args.index_regex, args.distance, args.distance_type
    )

    results = defaultdict(
        Counter, [(str(x), Counter()) for x in range(args.distance + 1)]
    )
    with dnaio.open(args.input) as reader:
        for record in reader:
            match = regex_pattern.search(
                record.sequence,
            )
            if match:
                fuzziness = str(sum(match.fuzzy_counts))
                results[fuzziness][(match.groups()[1:])] += 1

    if not results:
        logging.warning("No matches found.")
        return

    if args.output:
        write_results(results, args.output)
    else:
        for key, value in results.items():
            for seq, cnt in value.most_common():
                print(f"{key} {' '.join(list(seq))} {cnt}")
