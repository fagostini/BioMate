"""Blabber main script."""

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
    """Initialise module subparser."""
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
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--seq-length",
        type=int,
        default=100,
        help="Length of each sequence (default: 50)",
    )
    group.add_argument(
        "--seq-mask",
        type=str,
        default="Y100",
        help="Structure mask (default: Y100)",
    )
    parser.add_argument(
        "--index1",
        type=str,
        default=None,
        help="Index 1 sequence",
    )
    parser.add_argument(
        "--index2",
        type=str,
        default=None,
        help="Index 2 sequence",
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
        help="""Add to the undetermined files approximately 10%% of the sequences from each project
        in the same lane, if there are projects with different sequencing designs.
        """,
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="Seed for the random number generator (default: None, which uses system time)",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def parse_sequence_mask(
    mask_string: str, index1: str = None, index2: str = None
) -> dict:
    """
    Parse the sequence mask.
    """
    mask_split = re.findall(
        r"[^\W\d_]+|\d+", mask_string
    )  # Split characters from digits
    mask_dict = {"R1": 0, "I1": 0, "I2": 0, "R2": 0}  # Create empty dict
    key = "R1"  # Set the initial key value
    for i in range(len(mask_split)):
        if i % 2 == 0:  # Even positions should always be characters
            if mask_split[i] == "I":
                key = "I1" if mask_dict["I1"] == 0 else "I2"
            elif (
                mask_split[i] == "N"
            ):  # Ns should change the key only if it is set to I
                key = "R2" if key in ["I1", "I2"] else "R1"
            elif mask_split[i] == "Y":
                key = "R1" if mask_dict["I1"] == 0 and mask_dict["I2"] == 0 else "R2"
            else:
                raise ValueError(
                    f"Invalid character found in sequence mask: {mask_split[i]}"
                )
        else:  # Odd positions should always be digits
            if mask_split[i].isdigit():
                # If specified, use the index sequence instead of the length
                if index1 and key == "I1" and len(index1) == int(mask_split[i]):
                    mask_dict[key] = index1
                # If specified, use the index sequence instead of the length
                elif index2 and key == "I2" and len(index2) == int(mask_split[i]):
                    mask_dict[key] = index2
                else:
                    mask_dict[key] += int(mask_split[i])
            else:
                raise ValueError(
                    f"Sequence mask characters must be followed by digits! {mask_split[i]}"
                )
    return mask_dict


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if not args.alphabet.isalpha():
        raise argparse.ArgumentTypeError("Letters must be alphabetic characters.")
    if args.seq_length <= 0:
        raise argparse.ArgumentTypeError("Length must be a positive integer.")
    if args.seq_mask:
        _ = parse_sequence_mask(args.seq_mask.upper())
    if args.seq_number <= 0:
        raise argparse.ArgumentTypeError(
            "Number of sequences must be a positive integer."
        )
    if args.index1 and not args.index1.isalpha():
        raise argparse.ArgumentTypeError(
            "Index 1 must contain only alphabetic characters."
        )
    if args.index2 and not args.index2.isalpha():
        raise argparse.ArgumentTypeError(
            "Index 2 must contain only alphabetic characters."
        )
    if args.format == "fastq":
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


def generate_flowcell_id(instrument: str = "NovaSeq") -> str:
    """
    Generate the flowcell identifier.
    """
    if instrument == "NovaSeq":
        return "".join(
            [str(random.randint(10, 99))]
            + [
                string.ascii_letters[random.randint(len(string.ascii_letters))].upper()
                for i in range(6)
            ]
            + [str(random.randint(10))]
        )


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

    # Find the line where the 'BCLConvert_Data' or 'Data' section starts
    skip_lines = [
        i
        for i, line in enumerate(lines)
        if line.startswith("[BCLConvert_Data]") or line.startswith("[Data]")
    ]
    if not skip_lines:
        logging.error("No valid SampleSheet found! Data section not found in the file.")
        exit(1)
    else:
        skip_lines = skip_lines[0] + 1

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

    # Exit if the design is not present among the data columns
    if "OverrideCycles" not in data.columns and "Recipe" not in data.columns:
        raise RuntimeError("OverrideCycles or Recipe column not found in Data section!")

    # SampleSheet version 2
    if any([line.startswith("[BCLConvert_Data]") for line in lines]):
        # Create a dictionary to store the sample sheet data
        sample_sheet_dict = defaultdict(dict)
        counter = 1
        for proj in data.get_column("Sample_Project").unique(maintain_order=True):
            sample_sheet_dict.setdefault(proj, defaultdict(list))
            for lane, id, recipe, i1, i2 in (
                data.filter(polars.col("Sample_Project") == proj)
                .select("Lane", "Sample_ID", "OverrideCycles", "index", "index2")
                .iter_rows()
            ):
                sample_sheet_dict[proj].setdefault(
                    id,
                    tuple((f"S{counter}_{lane}", parse_sequence_mask(recipe, i1, i2))),
                )
                counter += 1
        return (generate_flowcell_id(), sample_sheet_dict)
    else:
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
                    id,
                    tuple(
                        (
                            f"S{counter}_{lane}",
                            parse_sequence_mask("Y" + recipe.replace("-", "Y")),
                        )
                    ),
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
    """Generate FASTQ sequence files using the dnaio library for I/O."""
    sampled_sequences = []
    read1_len = recipe["R1"]
    read2_len = recipe["R2"]
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
                    sequence="".join(random.choice(list(nucleotides), size=read1_len)),
                    qualities="I" * read1_len,
                )
            ]
            if read2_len != 0:
                reads.append(
                    dnaio.SequenceRecord(
                        name=f"{prefix}:{suffix}{extension.replace(' 1', ' 2')}",
                        sequence="".join(
                            random.choice(list(nucleotides), size=read2_len)
                        ),
                        qualities="I" * read2_len,
                    )
                )
            writer.write(*reads)
            if sample and random.randint(10) == 0:
                sampled_sequences.append(tuple(reads))
    return sampled_sequences


def main(args: argparse.Namespace) -> None:
    """Main function."""
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

        if sample_sheet:
            # unique_lanes = defaultdict(set)
            taint_sequences = defaultdict(list)
            for proj, samples in sample_sheet.items():
                for sample, (index, recipe) in samples.items():
                    lane_key = re.sub(r"S\d+_", "S0_", index)
                    # unique_lanes.setdefault(lane_key, set()).add(recipe)
                    index1 = recipe["I1"]
                    index2 = recipe["I2"]
                    if index1 != 0:
                        extension = f" 1:N:0:{index1}"
                        if index2 != 0:
                            extension += f"+{index2}"
                    else:
                        extension = (
                            f" 1:N:0:{'N' * 8}" if args.format == "fastq-ext" else ""
                        )
                    output_files = [
                        output_basepath.joinpath(proj)
                        .joinpath(sample)
                        .joinpath(
                            f"{sample.replace('Sample_', '')}_{index}_R{i}_001.fastq.gz"
                        )
                        for i, _ in enumerate(
                            {
                                k: v
                                for k, v in recipe.items()
                                if k in ["R1", "R2", "R3"]
                            },
                            start=1,
                        )
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
            # taint_sequences = {
            #     key: value
            #     for key, values in taint_sequences.items()
            #     if len(unique_lanes[key]) > 1
            #     for value in values
            #     if value
            # }
            # if args.taint and taint_sequences:
            #     for lane, sequences in taint_sequences.items():
            #         logging.debug(
            #             f"Adding {len(sequences)} sequences to the undetermined files for lane {lane}."
            #         )
            # for lane in unique_lanes:
            #     recipe = [
            #         int(x) for x in sorted(list(unique_lanes[lane]))[0].split("-")
            #     ]
            #     output_files = [
            #         output_basepath.joinpath(f"Undetermined_{lane}_R{i}_001.fastq.gz")
            #         for i, _ in enumerate(recipe, start=1)
            #     ]

            #     _ = generate_dnaio_fastq_files(
            #         output_files,
            #         nucleotides,
            #         args.seq_number,
            #         recipe,
            #         prefix,
            #         extension,
            #         sequences=taint_sequences.get(lane, []),
            #     )

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
