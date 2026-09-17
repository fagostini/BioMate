"""Blabber main script."""

import argparse
import contextlib
import logging
import pathlib
import re
import shutil
import string
import tempfile
from collections import defaultdict
from datetime import datetime
from itertools import product, batched
import dnaio
import polars
import numpy as np

from biomate.setup import setup_logging

TILE_WIDTH = 320  # 5120
TILE_HEIGHT = 160  # 2879
SURFACE_COUNT = 2
SWATH_COUNT = 4
TILE_COUNT = 16

# Module-level random number generator (initialized in main())
rng = np.random.default_rng()


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
        default=None,
        help="Structure mask (format: R1;I1;I2;R2)",
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
        "--flowcell-id",
        type=str,
        default=None,
        help="Use this flowcell id instead of generating a random string (default: None)",
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
    mask_string: str, index1: str | None = None, index2: str | None = None
) -> dict:
    """
    Parse the sequence mask.


    """
    mask_dict: dict[str, int | str] = {
        "U1": 0,
        "R1B": 0,
        "R1": 0,
        "R1A": 0,
        "I1B": 0,
        "I1": 0,
        "I1A": 0,
        "U2": 0,
        "I2B": 0,
        "I2": 0,
        "I2A": 0,
        "R2B": 0,
        "R2": 0,
        "R2A": 0,
    }  # Create empty dict

    mask_split = mask_string.split(";")

    if len(mask_split) > 4:
        raise ValueError(f"OverrideCycles mask has incorrect format: {mask_string}")

    for i_substring, substring in enumerate(mask_split):
        substring_type = "R" if "Y" in substring else "I"
        substring_index = i_substring // 2 + 1
        submask_split = re.findall(r"[^\W\d_]+|\d+", substring)

        for label, value in batched(submask_split, 2):
            suffix = ""
            if label == "U":
                substring_key = "U1" if mask_dict["U1"] == 0 else "U2"
                mask_dict[substring_key] = int(value)
            elif label in ["Y", "I"]:
                mask_dict[f"{substring_type}{substring_index}"] = int(value)
            else:
                suffix += (
                    "B" if mask_dict[f"{substring_type}{substring_index}"] == 0 else "A"
                )
                mask_dict[f"{substring_type}{substring_index}{suffix}"] = int(value)

    if index1:
        if len(index1) == mask_dict["I1"]:
            mask_dict["I1"] = index1
        else:
            raise ValueError(
                f"Error: Index 1 mask value ({mask_dict['I1']}) does not match the specified Index 1 length ({len(index1)})!"
            )

    if index2:
        if len(index2) == mask_dict["I2"]:
            mask_dict["I2"] = index2
        else:
            raise ValueError(
                f"Error: Index 2 mask value ({mask_dict['I2']}) does not match the specified Index 2 length ({len(index2)})!"
            )

    logging.debug(f"Mask '{mask_string}' parsed as {mask_dict}")
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
    return "".join(
        [str(rng.integers(10, 100))]
        + [
            string.ascii_letters[rng.integers(len(string.ascii_letters))].upper()
            for i in range(6)
        ]
        + [str(rng.integers(0, 10))]
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
        .with_columns(
            polars.col("Lane").cast(polars.String).str.zfill(3).alias("Lane_name")
        )
        .with_columns(polars.col("Lane_name").str.pad_start(4, "L"))
    )

    # Exit if the design is not present among the data columns
    if "OverrideCycles" not in data.columns and "Recipe" not in data.columns:
        raise RuntimeError("OverrideCycles or Recipe column not found in Data section!")

    # SampleSheet version 2
    if any([line.startswith("[BCLConvert_Data]") for line in lines]):
        # Create a dictionary to store the sample sheet data
        sample_sheet_dict = defaultdict(dict)
        counter = 1
        for proj in (
            data.get_column("Sample_Project").drop_nulls().unique(maintain_order=True)
        ):
            sample_sheet_dict.setdefault(proj, defaultdict(list))
            for lane, name, id, recipe, i1, i2 in (
                data.filter(polars.col("Sample_Project") == proj)
                .select(
                    "Lane",
                    "Lane_name",
                    "Sample_ID",
                    "OverrideCycles",
                    "index",
                    "index2",
                )
                .iter_rows()
            ):
                sample_sheet_dict[proj].setdefault(
                    id,
                    tuple(
                        (
                            f"S{counter}_{name}",
                            lane,
                            parse_sequence_mask(recipe, i1, i2),
                        )
                    ),
                )
                counter += 1
        return (generate_flowcell_id(), sample_sheet_dict)
    else:
        # Create a dictionary to store the sample sheet data
        sample_sheet_dict = defaultdict(dict)
        counter = 1
        for proj in data.get_column("Sample_Project").unique(maintain_order=True):
            sample_sheet_dict.setdefault(proj, defaultdict(list))
            for lane, name, id, recipe in (
                data.filter(polars.col("Sample_Project") == proj)
                .select("Lane", "Lane_name", "Sample_ID", "Recipe")
                .iter_rows()
            ):
                sample_sheet_dict[proj].setdefault(
                    id,
                    tuple(
                        (
                            f"S{counter}_{name}",
                            lane,
                            parse_sequence_mask("Y" + recipe.replace("-", "Y")),
                        )
                    ),
                )
                counter += 1
        return (data.get_column("FCID").unique().item(), sample_sheet_dict)


def generate_sequences_set(nucleotides: set, length: int, number: int) -> set:
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
    max_unique = len(nucleotides) ** length
    while len(sequences) < number and len(sequences) < max_unique:
        seq = "".join(rng.choice(list(nucleotides), size=length))
        sequences.add(seq)
    return sequences


def generate_dnaio_fastq_files(
    output_files: list[pathlib.Path],
    nucleotides: set,
    seq_number: int,
    recipe: dict[str, int | str],
    prefix: str,
    extension: str,
    tiles: list[str],
    taint: bool = False,
    sequences: list | None = None,
) -> list:
    """Generate FASTQ sequence files using the dnaio library for I/O.

    Args:
        output_files: List of output file paths
        nucleotides: Set of allowed nucleotides
        seq_number: Number of sequences to generate
        recipe: Recipe dict with R1, R2, U1, U2 lengths
        prefix: Read name prefix
        extension: Read name extension
        tiles: List of tile identifiers; sequences are assigned to them
            round-robin so that every tile is used at least once
        taint: If True, return ~10% of sequences for contamination simulation
        sequences: Pre-generated sequences to use

    Returns:
        List of sampled sequences for tainting (if taint=True), empty list otherwise
    """
    sampled_sequences = []
    k = int(np.ceil(np.sqrt(seq_number)))
    available_pos_x = [
        str(x)
        for x in rng.choice(range(1000, TILE_WIDTH * 10 + 1), size=k, replace=False)
    ]
    available_pos_y = [
        str(x)
        for x in rng.choice(range(1000, TILE_HEIGHT * 10 + 1), size=k, replace=False)
    ]
    assert tiles, "tiles must be provided"
    read1_len = int(recipe["R1"])
    read2_len = int(recipe["R2"])
    umi1 = int(recipe["U1"])
    umi2 = int(recipe["U2"])
    logging.debug(
        f"Generating sequence for {[x.name for x in output_files]} using tile(s): {tiles}"
    )
    with dnaio.open(*output_files, mode="w", fileformat="fastq") as writer:
        for reads in sequences or []:
            if reads:
                writer.write(*reads)
        for i in range(seq_number):
            suffix = ":".join(
                [
                    tiles[i % len(tiles)],
                    available_pos_x[i // k],
                    available_pos_y[i % k],
                ]
            )
            if umi1 != 0:
                suffix += ":" + "".join(rng.choice(list(nucleotides), size=umi1))
            if umi2 != 0:
                suffix += "+" if umi1 != 0 else ":"
                suffix += "".join(rng.choice(list(nucleotides), size=umi2))
            reads = [
                dnaio.SequenceRecord(
                    name=f"{prefix}:{suffix}{extension}",
                    sequence="".join(rng.choice(list(nucleotides), size=read1_len)),
                    qualities="I" * read1_len,
                )
            ]
            if read2_len != 0:
                reads.append(
                    dnaio.SequenceRecord(
                        name=f"{prefix}:{suffix}{extension.replace(' 1', ' 2')}",
                        sequence="".join(rng.choice(list(nucleotides), size=read2_len)),
                        qualities="I" * read2_len,
                    )
                )
            writer.write(*reads)

            # Collect ~10% of reads for taint contamination if enabled
            if taint and rng.random() < 0.1:
                sampled_sequences.append(reads)

    return sampled_sequences


def split(a, n):
    """Split a list into n approximately equal parts."""
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))


def partition_tiles_by_lane(available_tiles: list, lanes: list[str]) -> dict[str, list]:
    """Distribute the available tiles evenly across the lanes.

    Each lane receives a pool of tiles whose sizes differ by at most one.

    Args:
        available_tiles: List of tile identifiers.
        lanes: List of lane names.

    Returns:
        dict: Mapping of lane names to their tile pools.
    """
    if len(lanes) > len(available_tiles):
        raise ValueError(
            f"Cannot partition {len(available_tiles)} tiles across {len(lanes)} lanes: "
            "each lane must receive at least one tile."
        )
    return dict(zip(lanes, split(available_tiles, len(lanes))))


def tiles_per_lane(
    pool_size: int, lane_sample_counts: list[int], seq_number: int
) -> int:
    """Compute the number of tiles to use per lane.

    All tiles of a lane's pool are used unless some lane has fewer sequences
    than tiles, in which case the per-lane tile count is capped so that every
    lane uses the same number of tiles and every used tile gets at least one
    sequence.

    Args:
        pool_size: Number of tiles available to each lane.
        lane_sample_counts: Number of samples in each lane.
        seq_number: Number of sequences generated per sample.

    Returns:
        int: Number of tiles to use per lane.
    """
    return min(pool_size, min(count * seq_number for count in lane_sample_counts))


def assign_tiles_to_samples(tile_pool: list, n_samples: int) -> list[list]:
    """Distribute a lane's tile pool across its samples.

    Every sample receives at least one tile and every tile in the pool is
    assigned to at least one sample, so that the whole pool is used.

    Args:
        tile_pool: List of tile identifiers belonging to a lane.
        n_samples: Number of samples in the lane.

    Returns:
        list: List of tile lists, one per sample.
    """
    if n_samples <= 0:
        raise ValueError("n_samples must be a positive integer.")
    if not tile_pool:
        raise ValueError("tile_pool must not be empty.")
    pool_size = len(tile_pool)
    return [
        tile_pool[i::n_samples] or [tile_pool[i % pool_size]] for i in range(n_samples)
    ]


def generate_undetermined_files(
    output_basepath: pathlib.Path,
    taint_sequences: dict[str, list],
    taint_rate: float = 0.1,
) -> None:
    """Generate Undetermined FASTQ files with tainted sequences from other samples.

    Creates Undetermined_R1_001.fastq.gz and Undetermined_R2_001.fastq.gz in each lane
    directory, containing approximately taint_rate*100% reads from other samples.
    This simulates index hopping/cross-project contamination (~10% by default).

    Args:
        output_basepath: Base output path containing lane directories
        taint_sequences: Dict mapping lane_key to list of read samples from all projects
        taint_rate: Fraction of reads to include in undetermined files (default 0.1 = 10%)
    """
    for lane_key, sample_sequences in taint_sequences.items():
        # Extract lane number from key (format: S0_L001)
        lane_match = re.search(r"L(\d+)", lane_key)
        if not lane_match:
            logging.warning(f"Could not extract lane number from {lane_key}")
            continue

        lane_num = int(lane_match.group(1))

        # Pool the read pairs (lists of one or two SequenceRecords) from all samples
        # in this lane, keeping the R1/R2 pairing intact
        all_read_pairs = []
        for sample_list in sample_sequences:
            for reads in sample_list:
                all_read_pairs.append(reads if isinstance(reads, list) else [reads])

        if not all_read_pairs:
            logging.debug(f"No taint sequences for lane {lane_num}")
            continue

        # Shuffle and take taint_rate fraction
        rng.shuffle(all_read_pairs)
        n_taint = max(1, int(len(all_read_pairs) * taint_rate))
        taint_pool = all_read_pairs[:n_taint]

        logging.info(
            f"Writing {len(taint_pool)} tainted sequences to Undetermined files for lane {lane_num}"
        )

        # Write undetermined files to the parent directory
        undetermined_dir = output_basepath.parent / "Demultiplexing"
        undetermined_r1 = (
            undetermined_dir / f"Undetermined_S0_L{lane_num:03d}_R1_001.fastq.gz"
        )
        undetermined_r2 = (
            undetermined_dir / f"Undetermined_S0_L{lane_num:03d}_R2_001.fastq.gz"
        )

        undetermined_dir.mkdir(parents=True, exist_ok=True)

        # Write the R1 records to the R1 file and the R2 records (if any) to the
        # R2 file, using separate writers so that single-end reads are supported
        r2_present = any(len(reads) > 1 for reads in taint_pool)
        with dnaio.open(
            str(undetermined_r1), mode="w", fileformat="fastq"
        ) as writer_r1:
            r2_writer = (
                dnaio.open(str(undetermined_r2), mode="w", fileformat="fastq")
                if r2_present
                else contextlib.nullcontext()
            )
            with r2_writer as writer_r2:
                for reads in taint_pool:
                    writer_r1.write(reads[0])
                    if len(reads) > 1:
                        writer_r2.write(reads[1])


def main(args: argparse.Namespace) -> None:
    """Main function."""
    global rng
    if not args.output:
        args.quiet = False
    setup_logging(args)
    logging.info("Running blabber module...")

    rng = np.random.default_rng(args.random_seed)
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
                list(rng.choice([x for x in string.ascii_uppercase], size=2))
                + [str(rng.integers(10000, 100000))]
            )
            field2 = str(rng.integers(100, 1000))
            output_basepath = args.output.joinpath(
                f"{field0}_{field1}_{field2:>04}_A{field3}"
                if not args.flowcell_id
                else args.flowcell_id
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
            field1 = "".join(rng.choice(list(string.ascii_uppercase), size=2)) + str(
                rng.integers(10000, 100000)
            )
            field2 = str(rng.integers(100, 1000))
            field3 = str(rng.integers(10, 100)) + "".join(
                rng.choice(list(string.ascii_uppercase), size=6)
            )
            prefix = f"{field1}:{field2}:{field3}"
            extension = f" 1:N:0:{'N' * 8}" if args.format == "fastq-ext" else ""

        if sample_sheet:
            # unique_lanes = defaultdict(set)
            taint_sequences = defaultdict(list)
            available_tiles = [
                "".join(x)
                for x in product(
                    [str(i + 1) for i in range(SURFACE_COUNT)],
                    [str(i + 1) for i in range(SWATH_COUNT)],
                    [f"{i + 1:02}" for i in range(TILE_COUNT)],
                )
            ]
            rng.shuffle(available_tiles)

            lanes = sorted(
                {
                    lane
                    for samples in sample_sheet.values()
                    for (_, lane, _) in samples.values()
                }
            )
            lane_tile_pools = partition_tiles_by_lane(available_tiles, lanes)
            lane_samples = defaultdict(list)
            for proj, samples in sample_sheet.items():
                for sample, (_, lane, _) in samples.items():
                    lane_samples[lane].append((proj, sample))
            tiles_per_lane_count = tiles_per_lane(
                min(len(pool) for pool in lane_tile_pools.values()),
                [len(samples) for samples in lane_samples.values()],
                args.seq_number,
            )
            sample_tiles = {}
            for lane, tile_pool in lane_tile_pools.items():
                n_samples = len(lane_samples[lane])
                for (proj, sample), tiles in zip(
                    lane_samples[lane],
                    assign_tiles_to_samples(
                        tile_pool[:tiles_per_lane_count], n_samples
                    ),
                ):
                    sample_tiles[(proj, sample)] = tiles
            logging.info(
                f"Using {tiles_per_lane_count} tiles per lane "
                f"({tiles_per_lane_count * len(lanes)} of "
                f"{len(available_tiles)} available tiles)"
            )

            for proj, samples in sample_sheet.items():
                for sample, (index, lane, recipe) in samples.items():
                    proj_prefix = prefix + f":{lane}"  # Lane number
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
                            f"{sample.replace('Sample_', '')}_{index}_{read}_001.fastq.gz"
                        )
                        for read, length in recipe.items()
                        if read in ["R1", "R2", "R3"] and length != 0
                    ]
                    if not output_files[0].parent.is_dir():
                        output_files[0].parent.mkdir(parents=True, exist_ok=True)

                    sampled_sequences = generate_dnaio_fastq_files(
                        output_files,
                        nucleotides,
                        args.seq_number,
                        recipe,
                        proj_prefix,
                        extension,
                        sample_tiles[(proj, sample)],
                        args.taint,
                    )
                    taint_sequences.setdefault(lane_key, []).append(sampled_sequences)

            # Generate Undetermined FASTQ files with tainted sequences if --taint is enabled
            if args.taint and taint_sequences:
                generate_undetermined_files(output_basepath, taint_sequences)

        else:
            sequences = set()
            while len(sequences) < args.seq_number:
                seq = "".join(rng.choice(list(nucleotides), size=args.seq_length))
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
                            [str(rng.integers(0, 10000))]  # Tile number
                            + [":"]
                            + [str(rng.integers(0, 10000))]  # X coordinate
                            + [":"]
                            + [str(rng.integers(0, 10000))]  # Y coordinate
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
                        [str(rng.integers(0, 10000))]
                        + [":"]
                        + [str(rng.integers(0, 10000))]
                        + [":"]
                        + [str(rng.integers(0, 10000))]
                    )
                    print(
                        f"@{prefix}:{suffix}{extension}\n{seq}\n+\n{'I' * args.seq_length}"
                    )
    else:
        sequences = set()
        while len(sequences) < args.seq_number:
            seq = "".join(rng.choice(list(nucleotides), size=args.seq_length))
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
