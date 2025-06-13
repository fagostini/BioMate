import argparse
import logging
import pathlib
import re
import shutil
import string
import tempfile
from collections import defaultdict
from datetime import datetime

import dnaio
import polars
from numpy import random

from biomate.setup import setup_logging


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="FASTA, FASTQ or plain text sequence generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Generate random nucleotide sequences in FASTA, FASTQ or plain text format",
    )
    parser.add_argument(
        "--seq-number",
        type=int,
        default=100,
        help="Number of sequences to generate (default: 100)",
    )
    parser.add_argument(
        "--seq-length",
        type=int,
        default=100,
        help="Length of each sequence (default: 50)",
    )
    parser.add_argument(
        "--alphabet",
        type=str,
        default="ACGT",
        help="Alphabet of letters to use for generating sequences (default: 'ACGT')",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        default=None,
        help="Path to the output file where the sequences will be saved (default: stdout)",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["fasta", "fastq", "fastq-ext", "text"],
        default="text",
        help="Format of the output sequences (default: 'text')",
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        default=None,
        help="Path to a sample sheet file to use for generating sequences (default: None)",
    )
    parser.add_argument(
        "--taint",
        action="store_true",
        help="""Add to the undetermined files approximately 10%% of sequences from each project
        in the same lane, if there are projects with different sequencing designs.
        """,
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Seed for the random number generator (default: None, which uses system time)",
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
    if not args.alphabet.isalpha():
        raise argparse.ArgumentTypeError("Letters must be alphabetic characters.")
    if args.seq_length <= 0:
        raise argparse.ArgumentTypeError("Length must be a positive integer.")
    if args.seq_number <= 0:
        raise argparse.ArgumentTypeError(
            "Number of sequences must be a positive integer."
        )
    if "fastq" in args.format:
        if args.sample_sheet:
            if not args.sample_sheet.is_file():
                raise argparse.ArgumentTypeError(
                    f"Sample sheet file {args.sample_sheet} does not exist."
                )
            if not args.output:
                raise argparse.ArgumentTypeError(
                    "Output path is required when using a sample sheet."
                )
            else:
                if not args.output.is_dir():
                    raise argparse.ArgumentTypeError(
                        f"Output path {args.output} must be a directory."
                    )
    if args.output and not args.output.parent.is_dir():
        raise argparse.ArgumentTypeError(
            f"Output directory {args.output.parent} does not exist."
        )
    return args


def parse_sample_sheet(sample_sheet: pathlib.Path) -> tuple[str, dict]:
    """
    Parse the sample sheet file and return a dictionary of sample names and their corresponding lanes.

    Args:
        sample_sheet (pathlib.Path): Path to the sample sheet file.

    Returns:
        dict: Dictionary of sample names and their corresponding lanes.
    """
    # Read the sample sheet file and stores all lines in a list
    with open(sample_sheet) as input_file:
        lines = input_file.readlines()

    # Find the line where the FCID starts, where the actual data in CSV format starts
    skip_lines = [i for i, line in enumerate(lines) if line.startswith("FCID")]
    if not skip_lines:
        logging.error(
            "No valid SampleSheet found! 'FCID' header not found in the file. "
        )
        exit(1)
    else:
        skip_lines = skip_lines[0]

    # Write the lines to a temporary file, skipping the header lines
    tmp = tempfile.NamedTemporaryFile()
    with open(tmp.name, "w") as f:
        for line in lines[skip_lines:]:
            _ = f.write(line)

    # Read the temporary file using polars and process the data
    data = (
        polars.read_csv(tmp.name)
        .with_columns(polars.col("Lane").cast(polars.String).str.zfill(3))
        .with_columns(polars.col("Lane").str.pad_start(4, "L"))
    )

    # Create a dictionary to store the sample sheet data
    sample_sheet_dict = defaultdict(dict)
    counter = 1
    for proj in data.get_column("Sample_Project").unique(maintain_order=True):
        sample_sheet_dict.setdefault(proj, defaultdict(list))
        for lane, id, recipe in (
            data.filter(polars.col("Sample_Project") == proj)
            .select("Lane", "Sample_ID", "Recipe")
            .iter_rows()
        ):
            sample_sheet_dict[proj].setdefault(
                id, tuple((f"S{counter}_{lane}", recipe))
            )
            counter += 1

    return (data.get_column("FCID").unique().item(), sample_sheet_dict)


def generate_sequences_set(nucleotides: list, length: int, number: int) -> set:
    """
    Generate a set of random nucleotide sequences.

    Args:
        nucleotides (set): Set of nucleotides to use for generating sequences.
        length (int): Length of each sequence.
        number (int): Number of sequences to generate.

    Returns:
        set: A set of unique random nucleotide sequences.
    """
    sequences = set()
    while len(sequences) < number:
        seq = "".join(random.choice(list(nucleotides), size=length))
        sequences.add(seq)
    return sequences


def generate_dnaio_fastq_files(
    output_files: list[pathlib.Path],
    nucleotides: set,
    seq_number: int,
    recipe: list[int],
    prefix: str,
    extension: str,
    sample: bool = False,
    sequences: list = None,
) -> list:
    sampled_sequences = []
    with dnaio.open(*output_files, mode="w", fileformat="fastq") as writer:
        for reads in sequences or []:
            if reads:
                writer.write(*reads)
        for i in range(seq_number):
            suffix = "".join(
                [str(random.randint(10000))]
                + [":"]
                + [str(random.randint(10000))]
                + [":"]
                + [str(random.randint(10000))]
            )
            reads = [
                dnaio.SequenceRecord(
                    name=f"{prefix}:{suffix}{extension}",
                    sequence="".join(random.choice(list(nucleotides), size=recipe[0])),
                    qualities="I" * recipe[0],
                )
            ]
            if len(recipe) > 1:
                reads.append(
                    dnaio.SequenceRecord(
                        name=f"{prefix}:{suffix}{extension.replace(' 1', ' 2')}",
                        sequence="".join(
                            random.choice(list(nucleotides), size=recipe[1])
                        ),
                        qualities="I" * recipe[1],
                    )
                )
            writer.write(*reads)
            if sample and random.randint(10) == 0:
                sampled_sequences.append(tuple(reads))
    return sampled_sequences


def main(args: argparse.Namespace) -> None:
    if not args.output:
        args.quiet = False
    setup_logging(args)
    logging.info("Running blabber module...")

    random.seed(args.random_seed) if args.random_seed else random.seed()
    nucleotides = {x for x in args.alphabet.upper()}

    if args.format == "fasta":
        sequences = generate_sequences_set(
            nucleotides, args.seq_length, args.seq_number
        )
        if args.output:
            output_file = (
                args.output.joinpath("sequences.fasta")
                if args.output.is_dir()
                else args.output
            )
            with dnaio.open(output_file, mode="w", fileformat="fasta") as writer:
                for i, seq in enumerate(sequences, start=1):
                    writer.write(dnaio.SequenceRecord(name=f"seq{i}", sequence=seq))
        else:
            for i, seq in enumerate(sequences, start=1):
                print(f">seq{i}\n{seq}")

    elif args.format == "fastq" or args.format == "fastq-ext":
        if args.sample_sheet:
            field3, sample_sheet = parse_sample_sheet(args.sample_sheet)
            field0 = datetime.today().strftime("%Y%m%d")
            field1 = "".join(
                list(random.choice([x for x in string.ascii_uppercase], size=2))
                + [str(random.randint(10000, 99999))]
            )
            field2 = str(random.randint(100, 999))
            output_basepath = args.output.joinpath(
                f"{field0}_{field1}_{field2:>04}_A{field3}"
            )
            if not output_basepath.is_dir():
                output_basepath.mkdir(parents=True, exist_ok=True)

            shutil.copy(args.sample_sheet, output_basepath.joinpath("SampleSheet.csv"))

            output_basepath = output_basepath.joinpath("Demultiplexing")
            if not output_basepath.is_dir():
                output_basepath.mkdir(parents=True, exist_ok=True)

            prefix = f"{field1}:{field2}:{field3}"
        else:
            sample_sheet = None
            prefix = "".join(
                list(random.choice([x for x in string.ascii_uppercase], size=2))
                + [str(random.randint(10000, 99999))]
                + [":"]
                + [str(random.randint(100, 999))]
                + [":"]
                + [str(random.randint(10, 99))]
                + list(random.choice([x for x in string.ascii_uppercase], size=6))
            )

        prefix += "".join([":"] + [str(random.randint(10))])
        extension = f" 1:N:0:{'N' * 8}" if args.format == "fastq-ext" else ""

        if sample_sheet:
            unique_lanes = defaultdict(set)
            taint_sequences = defaultdict(list)
            for proj, samples in sample_sheet.items():
                for sample, (index, recipe) in samples.items():
                    lane_key = re.sub(r"S\d+_", "S0_", index)
                    unique_lanes.setdefault(lane_key, set()).add(recipe)
                    recipe = [int(x) for x in recipe.split("-")]
                    output_files = [
                        output_basepath.joinpath(proj)
                        .joinpath(sample)
                        .joinpath(
                            f"{sample.replace('Sample_', '')}_{index}_R{i}_001.fastq.gz"
                        )
                        for i, _ in enumerate(recipe, start=1)
                    ]
                    if not output_files[0].parent.is_dir():
                        output_files[0].parent.mkdir(parents=True, exist_ok=True)

                    sampled_sequences = generate_dnaio_fastq_files(
                        output_files,
                        nucleotides,
                        args.seq_number,
                        recipe,
                        prefix,
                        extension,
                        args.taint,
                    )
                    taint_sequences.setdefault(lane_key, []).append(sampled_sequences)
            taint_sequences = {
                key: value
                for key, values in taint_sequences.items()
                if len(unique_lanes[key]) > 1
                for value in values
                if value
            }
            if args.taint and taint_sequences:
                for lane, sequences in taint_sequences.items():
                    logging.debug(
                        f"Adding {len(sequences)} sequences to the undetermined files for lane {lane}."
                    )
            for lane in unique_lanes:
                recipe = [
                    int(x) for x in sorted(list(unique_lanes[lane]))[0].split("-")
                ]
                output_files = [
                    output_basepath.joinpath(f"Undetermined_{lane}_R{i}_001.fastq.gz")
                    for i, _ in enumerate(recipe, start=1)
                ]

                _ = generate_dnaio_fastq_files(
                    output_files,
                    nucleotides,
                    args.seq_number,
                    recipe,
                    prefix,
                    extension,
                    sequences=taint_sequences.get(lane, []),
                )

        else:
            sequences = set()
            while len(sequences) < args.seq_number:
                seq = "".join(random.choice(list(nucleotides), size=args.seq_length))
                sequences.add(seq)

            if args.output:
                output_file = (
                    args.output.joinpath("sequences.fastq")
                    if args.output.is_dir()
                    else args.output
                )
                with dnaio.open(output_file, mode="w", fileformat="fastq") as writer:
                    for i, seq in enumerate(sequences):
                        suffix = "".join(
                            [str(random.randint(10000))]
                            + [":"]
                            + [str(random.randint(10000))]
                            + [":"]
                            + [str(random.randint(10000))]
                        )
                        writer.write(
                            dnaio.SequenceRecord(
                                name=f"{prefix}:{suffix}{extension}",
                                sequence=seq,
                                qualities="I" * args.seq_length,
                            )
                        )
            else:
                for seq in sequences:
                    suffix = "".join(
                        [str(random.randint(10000))]
                        + [":"]
                        + [str(random.randint(10000))]
                        + [":"]
                        + [str(random.randint(10000))]
                    )
                    print(
                        f"@{prefix}:{suffix}{extension}\n{seq}\n+\n{'I' * args.seq_length}"
                    )
    else:
        sequences = set()
        while len(sequences) < args.seq_number:
            seq = "".join(random.choice(list(nucleotides), size=args.seq_length))
            sequences.add(seq)
        for seq in sequences:
            if args.output:
                output_file = (
                    args.output.joinpath("sequences.txt")
                    if args.output.is_dir()
                    else args.output
                )
                with open(output_file, "a") as writer:
                    writer.write(f"{seq}\n")
            else:
                print(seq)
