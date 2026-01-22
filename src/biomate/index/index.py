"""Index main script."""

import argparse
import logging
import pathlib
from collections import Counter, defaultdict
from itertools import product

import dnaio
import regex as re
from jellyfish import levenshtein_distance

from biomate.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
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
        default=0,
        help="Maximum distance between subject and query for sequences",
    )
    distance_group.add_argument(
        "--error-type",
        type=str,
        choices=["i", "d", "s", "e"],
        default="s",
        help="Type of error to consider: 'i' for insertion, 'd' for deletion, 's' for substitution, 'e' for any edit",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if not args.input.is_file():
        raise FileNotFoundError(f"Input file {args.input} does not exist.")
    return args


def compile_index_regex(index_regex: str, distance: int, error_type: str) -> re.Pattern:
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


def write_results(
    results: dict, total_records: int, output_path: pathlib.Path, file_index: int = 1
) -> None:
    """Write the results to the specified output path."""
    if not output_path.is_dir():
        output_path.mkdir(parents=True, exist_ok=True)

    max_len = max(
        [len(" ".join(list(seq))) for _, value in results.items() for seq in value]
    )

    with open(
        output_path.joinpath(f"pattern{file_index}_matches.txt"), "w"
    ) as output_file:
        for key, value in results.items():
            for seq, cnt in value.most_common():
                output_file.write(
                    f"{key} {' '.join(list(seq)):<{max_len}} {cnt} {cnt / total_records:.2%}\n"
                )

    with open(
        output_path.joinpath(f"pattern{file_index}_errors.txt"), "w"
    ) as output_file:
        for key, value in results.items():
            counts = sum(value.values())
            output_file.write(f"{key} {counts} {counts / total_records:.2%}\n")


def expand_regex(string: str) -> set[str]:
    """Expand a regex into all its string (genomic) patterns."""
    elements_list = []
    i = 0
    while i < len(string):
        c = string[i]
        if c == "." or c == "*":
            elements_list.append(["A", "C", "G", "T", "N"])
        elif c == "[":
            i += 1
            c = string[i]
            sublist = []
            while c != "]" and i < len(string):
                sublist.append(string[i])
                i += 1
                c = string[i]
            elements_list.append(sublist)
        elif c == "(":
            i += 1
            c = string[i]
            sublist = []
            substring = ""
            while c != ")" and i < len(string):
                if c == "|":
                    sublist.append(substring)
                    substring = ""
                else:
                    substring += string[i]
                i += 1
                c = string[i]
            sublist.append(substring)
            if string[i + 1] == "?":
                sublist.append("")
                i += 1
            elements_list.append(sublist)
        elif c == "+":
            raise ValueError(
                "The '+' character is currently not supported in the regex patterns is not supported!"
            )
        else:
            elements_list.append([string[i]])
        i += 1
    return set("".join(x) for x in product(*elements_list))


def validate_indexes_distance(patterns: list, distance: int) -> list[str, int]:
    """Check all indexes for compatible distances."""
    for query in patterns:
        expanded_query = expand_regex(query)
        for subject in patterns:
            print(f"Processing: {query} {subject}")
            if subject == query:
                continue
            expanded_subject = expand_regex(subject)
            if distance == 0:
                if any(
                    [
                        levenshtein_distance(q, s) <= distance
                        for q in expanded_query
                        for s in expanded_subject
                    ]
                ):
                    raise Exception(
                        f"Conflicting patterns '{query}' and '{subject}' at Levenshtein distance {distance}!"
                    )
                else:
                    print("No conflicts found!")
            else:
                all_ok = False
                valid_distance = distance
                while not all_ok and valid_distance >= 0:
                    all_ok = all(
                        [
                            levenshtein_distance(q, s) > valid_distance
                            for q in expanded_query
                            for s in expanded_subject
                        ]
                    )
                    print(
                        f"{expanded_query} {subject} {[levenshtein_distance(q, s) for q in expanded_query for s in expanded_subject]}"
                    )
                    if not all_ok:
                        valid_distance -= 1
                        print(f"Lowering distance to {valid_distance}")
                print(f"{all_ok} {valid_distance}")
                if valid_distance < 0:
                    raise Exception(
                        f"Conflicting patterns '{query}' and '{subject}' at Levenshtein distance 0!"
                    )
                else:
                    print(f"New distance: {valid_distance}")


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if not args.output:
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
    patterns_list = []
    if args.index_regex:
        logging.info(f"Using regex for indexing: {args.index_regex}")
        patterns_list.append(args.index_regex)
    if args.index_file:
        logging.info(f"Using index file: {args.index_file}")
        with open(args.index_file, "r") as index_file:
            for line in index_file:
                line = line.strip()
                if line:
                    patterns_list.append(line)

    for i, pattern in enumerate(patterns_list, start=1):
        logging.info(f"Processing pattern {i}: {pattern}")

        regex_pattern = compile_index_regex(pattern, args.distance, args.error_type)

        results = defaultdict(
            Counter, [(str(x), Counter()) for x in range(args.distance + 1)]
        )
        total_records = 0
        with dnaio.open(args.input) as reader:
            for record in reader:
                match = regex_pattern.search(
                    record.sequence,
                )
                if match:
                    fuzziness = str(sum(match.fuzzy_counts))
                    results[fuzziness][(match.groups()[1:])] += 1
                total_records += 1
        logging.debug(f"Total records processed: {total_records}")
        if not [x for x in results.values() if x]:
            logging.warning("   No matches found.")
            continue

        if args.output:
            write_results(results, total_records, args.output, file_index=i)
        else:
            for key, value in results.items():
                for seq, cnt in value.most_common():
                    print(
                        f"{key} {' '.join(list(seq))} {cnt} {cnt / total_records:.2%}"
                    )
