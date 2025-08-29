"""FastRewind main script."""

import argparse
import gzip
import logging
import pathlib
import struct
import tempfile
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime

import dnaio
import polars
import regex as re
from rich.progress import track

from biomate.setup import setup_logging

NOVASEQXPLUS_IMAGE_SIZE = (5120, 2879)


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description="Re-generate the BCL structure from demultiplexed flowcell folder",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="Re-generate the BCL structure from demultiplexed flowcell folder",
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        required=True,
        help="Path to the input flowcell directory containing demultiplexed FASTQ files",
    )
    parser.add_argument(
        "--output-path",
        type=pathlib.Path,
        required=False,
        default=pathlib.Path("./"),
        help="Path to the output directory where the BCL structure will be created (default: current directory)",
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        required=False,
        help="Path to the Sample Sheet file. If not provided, the script will look for SampleSheet.csv in the input directory.",
    )
    parser.add_argument(
        "--instrument",
        type=str,
        choices=["NovaSeqXPlus", "MiSeq", "NextSeq2000"],
        default="NovaSeqXPlus",
        help="Number of threads used by dnaio to read/write files. Default is 0, which corresponds to a single thread.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=0,
        help="Number of threads used by dnaio to read/write files. Default is 0, which corresponds to a single thread.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force overwrite of existing files in the output directory",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} is not a valid directory."
        )
    if not any(args.input_path.glob("**/*.fastq.gz")):
        raise argparse.ArgumentTypeError(
            f"No FASTQ files found in the input directory {args.input_path}."
        )
    if (
        not args.input_path.joinpath("SampleSheet.csv").is_file()
        and not args.sample_sheet
    ):
        raise argparse.ArgumentTypeError(
            f"SampleSheet.csv not found in the input directory {args.input_path}."
        )
    if args.sample_sheet and not args.sample_sheet.is_file():
        raise argparse.ArgumentTypeError(
            f"Sample sheet {args.sample_sheet} does not exist."
        )
    if args.output_path.is_dir() and not args.force:
        raise argparse.ArgumentTypeError(
            f"Output directory {args.output_path} already exists. Use --force to overwrite."
        )
    return args


def parse_flowcell_name(name: str) -> dict:
    """
    Parse the flowcell name to extract relevant information.
    """
    parts = name.split("_")
    if len(parts) < 4:
        raise ValueError(f"Invalid flowcell name format: {name}")
    return {
        "machine_id": parts[1],
        "run_id": parts[2],
        "flowcell_id": parts[3],
    }


def get_cycles(sample_sheet: pathlib.Path) -> "tuple[int, dict]":
    """
    Parse the sample sheet file and return the total number of cycles.

    Args:
        sample_sheet (pathlib.Path): Path to the sample sheet file.

    Returns:
        int: Total number of cycles calculated from the sample sheet.
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
    data = polars.read_csv(tmp.name).with_columns(
        [
            polars.col("Recipe").str.split_exact("-", 1).struct.unnest(),
            polars.col("index").str.len_chars().alias("index_length"),
            polars.col("index2").str.len_chars().alias("index2_length"),
        ]
    )
    if "field_1" not in data.columns:
        data = data.with_columns(polars.lit("0").alias("field_1"))
    data = data.with_columns(
        (
            polars.col("field_0").cast(polars.Int16)
            + polars.col("field_1").cast(polars.Int16)
            + polars.col("index_length").cast(polars.Int16)
            + polars.col("index2_length").cast(polars.Int16)
        ).alias("total_cycles")
    )
    cycles = (
        data.filter(polars.col("total_cycles") == polars.col("total_cycles").max())
        .head(1)
        .to_dicts()[0]
    )
    cycles_dict = {"R1": int(cycles["field_0"])}
    if "index_length" in cycles:
        cycles_dict.update({"I1": int(cycles["index_length"])})
    if "index2_length" in cycles:
        cycles_dict.update({"I2": int(cycles["index2_length"])})
    if "field_1" in cycles:
        cycles_dict.update({"R2": int(cycles["field_1"])})

    return (cycles["total_cycles"], cycles_dict)


def organise_fastq_files(path: pathlib.Path) -> OrderedDict:
    fastq_files = path.glob("**/*.fastq.gz")
    lanes_dict = defaultdict(defaultdict)
    for fastq_file in fastq_files:
        try:
            lane = re.search(r"_(L00[1-9])_", fastq_file.name).group(1)
            lanes_dict.setdefault(lane, defaultdict(list))
        except AttributeError:
            logging.error(
                f"Could not determine lane from file name {fastq_file.name}. "
                "Ensure the file name contains lane information in the format _LXXX_."
            )
            continue
        else:
            try:
                name = re.sub(r"_[RI][1-4]_", "_", fastq_file.name)
                name = name.replace(".fastq.gz", "")
            except Exception as e:
                logging.error(f"Error processing file {fastq_file.name}: {e}")
                continue
            else:
                lanes_dict[lane][name].append(fastq_file)
    return OrderedDict(sorted(lanes_dict.items()))


def parse_seqname_fields(name: str) -> dict:
    match = re.match(
        r"(?P<instrument>[A-Za-z0-9_]+):"
        r"(?P<run_number>[0-9]+):"
        r"(?P<flowcell_id>[-A-Za-z0-9]*):"
        r"(?P<lane>[0-9]+):"
        r"(?P<tile>[0-9]+):"
        r"(?P<x_pos>[0-9]+):"
        r"(?P<y_pos>[0-9]+)"
        r"( (?P<read>[12]):(?P<is_filtered>[YN]):(?P<control_number>[0-9]+):(?P<barcode_sequence>[A-Z]+(\+[A-Z]+)?))?",
        name,
    )
    return match.groupdict(match)


def generate_run_id(name: str) -> str:
    fields = parse_seqname_fields(name)
    run_id = (
        datetime.today().strftime("%y%m%d")
        + fields["instrument"]
        + "_"
        + fields["run_number"].zfill(4)
        + "_"
        + fields["flowcell_id"]
    )
    return run_id


def write_run_info_xml(
    xml_path: pathlib.Path,
    cycles_dict,
    lane_count,
    surface_count,
    swath_count,
    tile_count,
    tile_names,
    desc,
):
    xml_path.parent.mkdir(exist_ok=True, parents=True)
    fields = parse_seqname_fields(desc)
    run_id = generate_run_id(desc)
    run_number = fields["run_number"]
    flowcell_id = fields["flowcell_id"]
    instrument = fields["instrument"]
    date = datetime.now().strftime("%Y-%m-%dT%H:%M:%SZ")
    tile_count = int(tile_count / surface_count / swath_count)
    with open(xml_path, "wt") as f_out:
        f_out.write(
            '<?xml version="1.0"?> \n'
            '<RunInfo Version="6">\n'
            f'  <Run Id="{run_id}" Number="{run_number}">\n'
            f"    <Flowcell>{flowcell_id}</Flowcell>\n"
            f"    <Instrument>{instrument}</Instrument>\n"
            f"    <Date>{date}</Date>\n"
            f"    <Reads>\n"
        )
        for i, key in enumerate(cycles_dict.keys(), start=1):
            is_index = "N" if key in ["R1", "R2"] else "Y"
            f_out.write(
                f'      <Read Number="{i}" NumCycles="{cycles_dict[key]}" IsIndexedRead="{is_index}" IsReverseComplement="N"/>\n'
            )
        f_out.write(
            "    </Reads>\n"
            f'    <FlowcellLayout LaneCount="{lane_count}" SurfaceCount="{surface_count}" SwathCount="{swath_count}" TileCount="{tile_count}" />\n'
        )
        if tile_names:
            naming_convention = (
                "ThreeDigit" if len(tile_names[0].split("_")[1]) == 3 else "FourDigit"
            )
            f_out.write(
                f'      <TileSet TileNamingConvention="{naming_convention}">\n        <Tiles>\n'
            )
            for tile in tile_names:
                f_out.write(f"          <Tile>{tile}</Tile>\n")
            f_out.write("        <Tiles>\n      </TileSet>\n")
        f_out.write("    </FlowcellLayout>\n")
        f_out.write(
            f'    <ImageDimensions Width="{NOVASEQXPLUS_IMAGE_SIZE[0]}" Height="{NOVASEQXPLUS_IMAGE_SIZE[1]}"/>\n'
        )
        f_out.write(
            "    <ImageChannels>\n      <Name>blue</Name>\n      <Name>green</Name>\n    </ImageChannels>\n"
        )

        f_out.write("  </Run>\n</RunInfo>\n")


def write_filter(filter_path: pathlib.Path, cluster_count: int):
    filter_path.parent.mkdir(exist_ok=True, parents=True)
    with open(filter_path, "wb") as f_out:
        f_out.write(struct.pack("<III", 0, 3, cluster_count))
        f_out.write(bytes([1] * cluster_count))  # bit 0 is pass or failed filter


def write_control(control_path: pathlib.Path, cluster_count: int):
    control_path.parent.mkdir(exist_ok=True, parents=True)
    with open(control_path, "wb") as f_out:
        f_out.write(struct.pack("<III", 0, 3, cluster_count))
        # two bytes for each cluster
        f_out.write(bytes([0, 0] * cluster_count))


def encode_loc_bytes(x_pos, y_pos):
    # In the read name, the smallest value for X and Y is 1000, so they are
    # re-scaled in order for them to start from 0, while the division by 10 for
    # converting them to float (not sure if needed or if the factor is correct)
    return struct.pack("<ff", (int(x_pos) - 1000) / 10, (int(y_pos) - 1000) / 10)


def write_locs(locs_path: pathlib.Path, unique_locs: set):
    locs_path.parent.mkdir(exist_ok=True, parents=True)
    # The locs file format is one 3 Illumina formats(pos, locs, and clocs) that stores position data exclusively.
    # locs files store position data for successive clusters in 4 byte float pairs, described as follows:
    #     bytes 1-4    : (int?) Version number (1)
    #     bytes 5-8    : 4 byte float equaling 1.0
    #     bytes 9-12   : unsigned int numClusters
    #     bytes 13-16: : X coordinate of first cluster (little-endian 32-bit float)
    #     bytes 17-20: : Y coordinate of first cluster (little-endian 32-bit float)
    #     The remaining bytes of the file store the X and Y coordinates of the remaining clusters.
    with open(locs_path, "wb") as f_out:
        f_out.write(struct.pack("<IfI", 1, 1, len(unique_locs)))
        for x_pos, y_pos in sorted(unique_locs, key=lambda x: (x[1], x[0])):
            f_out.write(encode_loc_bytes(x_pos, y_pos))


def extract_flowcell_layout(lane_tiles: dict) -> tuple:
    df = (
        polars.from_dict({k: list(v) for k, v in lane_tiles.items()})
        .unpivot()
        .with_columns(
            [
                polars.col("value").cast(polars.UInt16),
                polars.lit("A").alias("Flowcell"),
                polars.col("variable").str.slice(3, 1).cast(polars.Int8).alias("Lane"),
                polars.col("value").str.slice(0, 1).cast(polars.Int8).alias("Surface"),
                polars.col("value").str.slice(1, 1).cast(polars.Int8).alias("Swath"),
            ]
        )
        .sort(["Lane", "value"])
        .with_columns(
            polars.concat_str(
                polars.col("Lane"), polars.col("value"), separator="_"
            ).alias("Name")
        )
    )

    return (
        tuple(
            (
                df.group_by("Flowcell")
                .agg(
                    [
                        polars.col("Lane").n_unique(),
                        polars.col("Surface").n_unique(),
                        polars.col("Swath").n_unique(),
                        polars.len().alias("Tile"),
                    ]
                )
                .drop("Flowcell")
                .transpose()
                .get_column("column_0")
                .to_list()
            )
        ),
        df.get_column("Name").unique(maintain_order=True).to_list(),
    )


def pack_cbcl_header(cycle: int, tiles_metrics: dict, quality_map: dict):
    h_version = 1
    # h_size = 6321 # Dynamically extracted based on the available tiles
    h_basebits = 2
    h_qbits = 2
    h_bins = 4
    return (
        struct.pack(
            "<HIBBI",
            h_version,
            48 + 16 * len(tiles_metrics) + 1,
            h_basebits,
            h_qbits,
            h_bins,
        )
        + struct.pack(
            "<" + "I" * len(quality_map.keys()) * 2,
            *[x for y in enumerate(quality_map.keys()) for x in y],
        )
        + struct.pack("<I", len(tiles_metrics))
        + b"".join(
            [
                struct.pack(
                    "<IIII",
                    int(tile_id),
                    tiles_metrics[tile_id]["clusters"],
                    tiles_metrics[tile_id]["uncompressed"][cycle],
                    tiles_metrics[tile_id]["compressed"][cycle],
                )
                for tile_id in tiles_metrics.keys()
            ]
        )
        + struct.pack("<B", 0)
    )


def encode_cluster_bits(
    base: str,
    qual: str,
    sequence_map: dict = {"A": "00", "C": "01", "G": "10", "T": "11"},
    quality_map: dict = {2: "00", 9: "01", 24: "10", 40: "11"},
    phred_offset: int = 33,
):
    return sequence_map.get(base, "00") + quality_map.get(
        ord(qual) - phred_offset,
        "00",
    )


def write_bcls(
    output_path: pathlib.Path,
    lane: str,
    tiles_path: pathlib.Path,
    total_cycles: int,
    threads: int = 0,
    instrument: str = "NovaSeqXPlus",
    quality_map={2: "00", 9: "01", 24: "10", 40: "11"},
):
    pattern = (
        f"{lane}/[1-2]_[1-2][1-2][0-9][0-9].fastq.gz"
        if instrument == "NovaSeqXPlus"
        else f"{lane[-1]}_*.fastq.gz"
    )

    # Get the list of all tile files, and group them by tile id
    tiles_files = defaultdict(list)
    for fq in tiles_path.glob(pattern):
        tiles_files.setdefault(fq.name.split("_")[0], list()).append(fq)

    # Cycle through surfaces and associated tiles
    for surface, fq_list in tiles_files.items():
        # Create the final output files list

        logging.debug(f"Processing Lane {lane} surface {surface}:")
        # Sort the files by tile number
        fq_sorted = sorted(
            fq_list,
            key=lambda fq: int(fq.name.split("_")[1].replace(".fastq.gz", "")),
        )
        # Create a dict with tile metrics (needed for the final file)
        tiles_metrics = defaultdict(dict)
        # Initialise a temporary directory where to put all the intermediate files
        # This will be deleted every time the processing of a surface is complete
        with tempfile.TemporaryDirectory(delete=True) as tempdir:
            # Create a dictionary with cycle files for each tile
            temp_files = {
                tile: [
                    pathlib.Path(tempdir).joinpath(f"{tile}_{cycle}.bcl")
                    for cycle in range(1, total_cycles + 1)
                ]
                for tile in [
                    fq.name.split("_")[1].replace(".fastq.gz", "") for fq in fq_sorted
                ]
            }
            # Cycle through the tile files
            for fq in fq_sorted:
                logging.debug(f"   {fq.name}")
                try:
                    tiles_size = 0
                    # Get the tile id
                    tile_id = fq.name.split("_")[1].replace(".fastq.gz", "")
                    tiles_metrics.setdefault(tile_id, dict()).update(
                        {"clusters": 0, "uncompressed": {}, "compressed": {}}
                    )
                    # Open all the cycle files for the given tile
                    opened_files = [open(file, "bw") for file in temp_files[tile_id]]
                    with dnaio.open(fq, open_threads=threads) as reader:
                        # Initialise a buffer for the bits. This is necessary because every read accounts
                        # only for half byte, so we need at least two to create a byte. One can probably read
                        # two at the time, but I thought that that would be more messy
                        buffer = []
                        for i, read in enumerate(reader, start=1):
                            # If the byte is complete
                            if i % 2 == 0:
                                # Cycle through the sequence and quality positions and extract the corresponding bits,
                                # adding them to those already stored in the buffer
                                for cycle in range(total_cycles):
                                    bits_string = buffer[cycle] + encode_cluster_bits(
                                        read.sequence[cycle], read.qualities[cycle]
                                    )
                                    # Convert the bits to a byte and write them to file
                                    opened_files[cycle].write(
                                        bytes([int(bits_string, 2)])
                                    )
                                buffer = []
                            else:
                                # Cycle through the sequence and quality positions and extract the corresponding bits
                                for cycle in range(total_cycles):
                                    buffer.append(
                                        encode_cluster_bits(
                                            read.sequence[cycle], read.qualities[cycle]
                                        )
                                    )
                            tiles_size += 1
                        # Add the number of clusters to the metrics
                        tiles_metrics[tile_id]["clusters"] = tiles_size
                        # If there is an odd number of sequences, the buffer would be incomplete
                        # Fill the rest of the buffer and write it to file
                        if buffer:
                            logging.debug(f"Left-over buffer: {buffer}")
                            for cycle in range(total_cycles):
                                opened_files[cycle].write(
                                    bytes([int(buffer[cycle].ljust(8, "0"), 2)])
                                )

                except RuntimeError as e:
                    logging.exception(
                        f"Warning: something went wrong! Quality string: {read.qualities} ({read.qualities[cycle]})"
                    )
                    logging.exception(e)
                finally:
                    logging.debug("Closing opened files...")
                    # Close all the opened files
                    for file in opened_files:
                        try:
                            file.close()
                        except IOError:
                            pass

            for cycle in range(total_cycles):
                with open(
                    pathlib.Path(tempdir).joinpath(f"all_tiles_{cycle + 1}"), "wb"
                ) as out_handle:
                    for tile_id in temp_files.keys():
                        if temp_files[tile_id][cycle].exists():
                            with open(temp_files[tile_id][cycle], "rb") as in_handle:
                                buffer = in_handle.read()
                                tiles_metrics[tile_id]["uncompressed"].update(
                                    {cycle: len(buffer)}
                                )
                                # Compress the lane bytes
                                buffer = gzip.compress(buffer, compresslevel=5)
                                tiles_metrics[tile_id]["compressed"].update(
                                    {cycle: len(buffer)}
                                )
                                out_handle.write(buffer)

            # for k, v in tiles_metrics.items():
            #     for cycle in range(total_cycles):
            #         print(
            #             k,
            #             cycle,
            #             v["clusters"],
            #             v["uncompressed"][cycle],
            #             v["compressed"][cycle],
            #         )

            for cycle in range(total_cycles):
                # Generate the header for the file
                header = pack_cbcl_header(cycle, tiles_metrics, quality_map)
                # Set and create the cycle folders
                cycle_path = output_path.joinpath(
                    f"Data/Intensities/BaseCalls/{lane}/C{cycle + 1}.1"
                )
                cycle_path.mkdir(exist_ok=True, parents=True)
                # Write header and body to the final file
                with open(
                    cycle_path.joinpath(f"{lane}_{surface}.cbcl"), "wb"
                ) as out_handle:
                    with open(
                        pathlib.Path(tempdir).joinpath(f"all_tiles_{cycle + 1}"), "rb"
                    ) as in_handle:
                        out_handle.write(header)
                        out_handle.write(in_handle.read())


def parse_fastq_groups(
    name: str,
    filenames: list,
    tempdir: pathlib.Path,
    threads: int = 0,
    instrument: str = "NovaSeqXPlus",
) -> dict:
    # Tiles list, used to update the Counter
    buffer_size = 1000000
    buffer_container = defaultdict(list)
    buffer_counter = 0
    tiles_and_locs = defaultdict(set)

    def process_read_name(name: str) -> tuple[str, list, str]:
        split_name = name.split(" ")
        indexes = (
            split_name[1].split(":")[-1].replace("+", "") if len(split_name) > 1 else ""
        )
        short_name = split_name[0]
        split_name = split_name[0].split(":")
        return (short_name, split_name, indexes)

    def dump_and_empty_buffer(
        buffer_container, output_path: pathlib.Path, threads: int = 0
    ) -> defaultdict:
        for key, values in buffer_container.items():
            write_mode = (
                "a" if output_path.joinpath(f"{key}.fastq.gz").exists() else "w"
            )
            with dnaio.open(
                output_path.joinpath(f"{key}.fastq.gz"),
                mode=write_mode,
                open_threads=threads,
            ) as writer:
                for read in values:
                    writer.write(read)
        return defaultdict(list)

    if instrument == "NovaSeqXPlus":
        tempdir = tempdir.joinpath(name)
        tempdir.mkdir(parents=True, exist_ok=True)

    with dnaio.open(*sorted(filenames), open_threads=threads) as reader:
        if len(filenames) == 1:
            for r1 in track(reader):
                r1.name, split_name, indexes = process_read_name(r1.name)
                # Assume the index is at the end
                r1.sequence += indexes
                r1.qualities += "I" * len(indexes)

                formatted_key = f"{split_name[3]}_{split_name[4]}"

                buffer_container.setdefault(formatted_key, list()).append(r1)
                buffer_counter += 1
                if buffer_counter % buffer_size == 0:
                    buffer_container = dump_and_empty_buffer(
                        buffer_container, tempdir, threads
                    )

                tiles_and_locs.setdefault(split_name[4], set()).add(
                    (split_name[5], split_name[6])
                )
        elif len(filenames) == 2:
            for r1, r2 in track(reader):
                r1.name, split_name, indexes = process_read_name(r1.name)
                r1.sequence += indexes + r2.sequence
                r1.qualities += "I" * len(indexes) + r2.qualities

                formatted_key = f"{split_name[3]}_{split_name[4]}"

                buffer_container.setdefault(formatted_key, list()).append(r1)
                buffer_counter += 1
                if buffer_counter % buffer_size == 0:
                    buffer_container = dump_and_empty_buffer(
                        buffer_container, tempdir, threads
                    )

                tiles_and_locs.setdefault(split_name[4], set()).add(
                    (split_name[5], split_name[6])
                )
        # For cases where more than 2 files are present, assume that the extra files are the indexes
        elif len(filenames) == 3:
            for r1, r2, r3 in track(reader):
                r1.name, split_name, indexes = process_read_name(r1.name)
                r1.sequence += r2.sequence + r3.sequence
                r1.qualities += r2.qualities + r3.qualities

                formatted_key = f"{split_name[3]}_{split_name[4]}"

                buffer_container.setdefault(formatted_key, list()).append(r1)
                buffer_counter += 1
                if buffer_counter % buffer_size == 0:
                    buffer_container = dump_and_empty_buffer(
                        buffer_container, tempdir, threads
                    )

                tiles_and_locs.setdefault(split_name[4], set()).add(
                    (split_name[5], split_name[6])
                )
        else:
            for r1, r2, r3, r4 in track(reader):
                r1.name, split_name, indexes = process_read_name(r1.name)
                r1.sequence += r2.sequence + r3.sequence + r4.sequence
                r1.qualities += r2.qualities + r3.qualities + r4.qualities

                formatted_key = f"{split_name[3]}_{split_name[4]}"

                buffer_container.setdefault(formatted_key, list()).append(r1)
                buffer_counter += 1
                if buffer_counter % buffer_size == 0:
                    buffer_container = dump_and_empty_buffer(
                        buffer_container, tempdir, threads
                    )

                tiles_and_locs.setdefault(split_name[4], set()).add(
                    (split_name[5], split_name[6])
                )

    if buffer_container:
        for key, values in buffer_container.items():
            write_mode = "a" if tempdir.joinpath(f"{key}.fastq.gz").exists() else "w"
            with dnaio.open(
                tempdir.joinpath(f"{key}.fastq.gz"),
                mode=write_mode,
                open_threads=threads,
            ) as writer:
                for read in values:
                    writer.write(read)

    for d in tempdir.iterdir():
        print(d)

    logging.info(f"Correctly processed {buffer_counter} clusters from sample '{name}'")
    return tiles_and_locs


def main(args: argparse.Namespace) -> None:
    """Main function."""
    setup_logging(args)
    logging.info("Running fastrewind module...")

    # Get individual and total number of cycles
    total_cycles, cycles_dict = get_cycles(
        args.sample_sheet or args.input_path.joinpath("SampleSheet.csv")
    )

    logging.debug("Gathering FASTQ files...")
    # Search and split the FASTQ files by lane and sample name
    lanes_dict = organise_fastq_files(args.input_path)

    # Get one read name to extract the run details
    sample_read = ""
    sample_lane = list(lanes_dict.keys())[0]
    sample_name = list(lanes_dict[sample_lane].keys())[0]
    sample_file = sorted(lanes_dict[sample_lane][sample_name])[0]
    with dnaio.open(sample_file, open_threads=args.threads) as reader:
        for read in track(reader):
            sample_read = read.name
            break

    # The next line is for debugging purposes
    outer_break = False
    # Dictionary containing all unique tiles per lane
    tiles_by_lane = defaultdict(set)
    flowcell_unique_locs = set()
    for lane, files in lanes_dict.items():
        logging.info(f"Processing lane '{lane}'")
        # Create a temporary directory where to put all intermediate files
        with tempfile.TemporaryDirectory(delete=True) as tempdir:
            tempdir = pathlib.Path(tempdir)
            # Dictionary containing the cluster counts by tile
            unique_locs_by_tile = defaultdict(set)
            for name, filenames in files.items():
                logging.info(f"   Processing '{name}': {[x.name for x in filenames]}")
                logging.debug("Parsing FASTQ files...")
                tiles_and_locs = parse_fastq_groups(
                    lane, filenames, tempdir, args.threads, args.instrument
                )
                # tiles_and_locs = parse_fastq_groups(
                #     lane, filenames, args.output_path, args.threads, args.instrument
                # )

                for k, v in tiles_and_locs.items():
                    unique_locs_by_tile.setdefault(k, set()).update(v)
                    flowcell_unique_locs.update(v)

                # Debugging code
                outer_break = True
                break

            logging.debug("Writing (C)BCL files...")
            write_bcls(
                args.output_path,
                lane,
                tempdir,
                total_cycles,
                args.threads,
                args.instrument,
            )

            # lane_file = pathlib.Path(tempdir).joinpath(f"{lane}.fastq.gz")
            # lane_file = args.output_path.joinpath(f"{lane}.fastq.gz")

            # try:
            #     lane_file.unlink()
            # except FileNotFoundError as e:
            #     logging.exception(f"Failed to remove file '{lane_file}'")

            tiles_by_lane.setdefault(lane, set()).update(
                {x for x in tiles_and_locs.keys()}
            )

            logging.debug("Writing FILTER files...")
            for tile, cc in unique_locs_by_tile.items():
                # Naming convention: Data/Intensities/BaseCalls/<L00_prefixed_lane>/s_<single_digit_lane>_<tile>.filter
                write_filter(
                    args.output_path.joinpath(
                        f"Data/Intensities/BaseCalls/{lane}/s_{lane[-1]}_{tile}.filter"
                    ),
                    cluster_count=len(cc),
                )
                # # Naming convention: Data/Intensities/BaseCalls/<L00_prefixed_lane>/s_<single_digit_lane>_<tile>.control
                # write_control(
                #     args.output_path.joinpath(
                #         f"Data/Intensities/BaseCalls/{lane}/s_{lane[-1]}_{tile}.control"
                #     ),
                #     cluster_count=cc,
                # )

            if outer_break:
                break

            logging.debug("Writing LOCS file(s)...")
        # Naming convention: Data/Intensities/BaseCalls/<L00_prefixed_lane>/s_<single_digit_lane>_<tile>.locs
        write_locs(
            args.output_path.joinpath("Data/Intensities/s.locs"),
            flowcell_unique_locs,
        )

    (lane_count, surface_count, swath_count, tile_count), tile_names = (
        extract_flowcell_layout(tiles_by_lane)
    )

    write_run_info_xml(
        args.output_path.joinpath("RunInfo.xml"),
        cycles_dict,
        lane_count,
        surface_count,
        swath_count,
        tile_count,
        tile_names,
        sample_read,
    )
