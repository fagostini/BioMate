"""Index main script."""

import argparse
import logging
import pathlib
import sys
from collections import Counter, defaultdict
from itertools import combinations, product
from math import comb
from multiprocessing import Pool

import dnaio
import regex as re

from biomate.headers import count_header_indexes
from biomate.setup import setup_logging

try:
    import biomate_kernels
except ImportError:  # optional acceleration kernel
    biomate_kernels = None

# Number of sequences passed to the rust scan engine per call
_KERNEL_CHUNK = 100_000


def _detect_input_type(filename: str) -> str | None:
    """Detect the input file type based on filename."""
    if re.match(r".*\.f(ast)?q(\.gz)?$", filename):
        return "FASTQ"
    elif re.match(r".*\.f(ast)?a(.gz)?$", filename):
        return "FASTA"
    return None


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
        nargs="+",
        required=True,
        help="Path(s) to the input FASTA/FASTQ file(s); results are kept per file",
    )
    parser.add_argument(
        "--output",
        type=pathlib.Path,
        required=False,
        help="Path to the output path where the results will be saved",
    )
    index_group = parser.add_mutually_exclusive_group(required=False)
    index_group.add_argument(
        "--index-regex",
        type=str,
        default=None,
        help="Regular expression to search in the FastQ file. Without it and "
        "--index-file, the indexes in the read headers are counted and grouped",
    )
    index_group.add_argument(
        "--index-file",
        type=pathlib.Path,
        default=None,
        help="Path to the index file. One index per line. Without it and "
        "--index-regex, the indexes in the read headers are counted and grouped",
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
    parser.add_argument(
        "--engine",
        type=str,
        choices=["regex", "rust"],
        default="regex",
        help="Search engine: 'regex' (default) or 'rust' (requires the optional biomate_kernels package)",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=1,
        help="Index discovery only: ignore header indexes seen fewer times "
        "(their total is reported separately). Default: 1, keep all",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of input files processed in parallel (with --output only)",
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    Validate the command line arguments.
    """
    if args.min_count < 1:
        raise argparse.ArgumentTypeError("--min-count must be >= 1.")
    if args.threads < 1:
        raise argparse.ArgumentTypeError("--threads must be >= 1.")
    for path in args.input:
        if not path.is_file():
            raise argparse.ArgumentTypeError(f"Input file {path} does not exist.")
    names = [_output_name(path) for path in args.input]
    if len(set(names)) != len(names):
        raise argparse.ArgumentTypeError(
            "Input files must have distinct names: per-file results would overwrite each other."
        )
    return args


def _output_name(path: pathlib.Path) -> str:
    """Input file name without the .gz and FASTA/FASTQ suffixes."""
    return re.sub(r"\.f(ast)?[aq](\.gz)?$", "", path.name)


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


def _sorted_by_count(counter: Counter) -> list[tuple[tuple[str, ...], int]]:
    """Counter items by descending count, then ascending group text.

    ``most_common()`` keeps first-seen order among ties, which differs between
    the regex path (insertion order) and the rust path (hash-map iteration),
    so the tie order would not be reproducible across engines or runs.
    """
    return sorted(counter.items(), key=lambda kv: (-kv[1], " ".join(list(kv[0]))))


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
            for seq, cnt in _sorted_by_count(value):
                output_file.write(
                    f"{key} {' '.join(list(seq)):<{max_len}} {cnt} {cnt / total_records:.2%}\n"
                )

    with open(
        output_path.joinpath(f"pattern{file_index}_errors.txt"), "w"
    ) as output_file:
        for key, value in results.items():
            counts = sum(value.values())
            output_file.write(f"{key} {counts} {counts / total_records:.2%}\n")


_BASES = "ACGTN"


def _neighbours(index: str, distance: int):
    """Yield ``(variant, k)`` for every string ``k`` (1..distance) substitutions
    away from ``index``, using A/C/G/T/N and keeping ``+`` separators fixed."""
    positions = [i for i, c in enumerate(index) if c != "+"]
    for k in range(1, distance + 1):
        for pos in combinations(positions, k):
            for subs in product(_BASES, repeat=k):
                if any(index[p] == b for p, b in zip(pos, subs)):
                    continue
                chars = list(index)
                for p, b in zip(pos, subs):
                    chars[p] = b
                yield "".join(chars), k


def group_indexes(
    counts: Counter, distance: int
) -> dict[str, list[tuple[int, str, int]]]:
    """Greedily group indexes within ``distance`` substitutions.

    The most frequent index seeds a group; each less frequent index joins the
    first (most frequent) seed of the same length within the Hamming distance,
    otherwise it seeds a new group. Members are ``(distance, index, count)``.
    """
    groups: dict[str, list[tuple[int, str, int]]] = {}
    rank: dict[str, int] = {}
    for index, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        best = None
        if distance > 0:
            n_variants = sum(
                comb(len(index), k) * (len(_BASES) - 1) ** k
                for k in range(1, distance + 1)
            )
            # Each index costs min(seeds, variants): scanning the seeds is
            # cheap early on, looking up the variants once seeds pile up
            if len(groups) <= n_variants:
                for seed in groups:
                    if len(seed) == len(index):
                        dist = sum(a != b for a, b in zip(seed, index))
                        if dist <= distance:
                            best = (dist, seed)
                            break
            else:
                hits = [
                    (rank[v], k, v)
                    for v, k in _neighbours(index, distance)
                    if v in rank
                ]
                if hits:
                    _, dist, seed = min(hits)
                    best = (dist, seed)
        if best is None:
            rank[index] = len(rank)
            groups[index] = [(0, index, count)]
        else:
            groups[best[1]].append((best[0], index, count))
    return groups


def _run_discovery(
    path: pathlib.Path, args: argparse.Namespace, out_dir: pathlib.Path | None
) -> None:
    """Count and group the header indexes of one input file."""
    if args.error_type != "s":
        logging.warning("Index discovery groups by substitutions only")
    try:
        counts = count_header_indexes(path)
    except ValueError as e:
        logging.warning(e)
        counts = Counter()
    total_records = sum(counts.values())
    if not counts:
        logging.warning("   No indexes found.")
        return
    # Rare indexes (mostly sequencing noise) are dropped before grouping:
    # it bounds grouping time, memory and output size on large files
    below_reads = below_indexes = 0
    if args.min_count > 1:
        kept = Counter()
        for index, count in counts.items():
            if count >= args.min_count:
                kept[index] = count
            else:
                below_reads += count
                below_indexes += 1
        counts = kept
    groups = group_indexes(counts, args.distance)
    del counts
    ordered = sorted(
        groups.items(), key=lambda kv: (-sum(c for _, _, c in kv[1]), kv[0])
    )
    below = (
        f"{below_reads} reads in {below_indexes} indexes below --min-count "
        f"{args.min_count} ({below_reads / total_records:.2%})"
    )
    if below_indexes:
        logging.info(below)

    if out_dir is None:
        _write_discovery(ordered, total_records, sys.stdout)
        if below_indexes:
            print(f"# {below}")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir.joinpath("indexes_matches.txt"), "w") as matches:
        _write_discovery(ordered, total_records, matches)
    with open(out_dir.joinpath("indexes_groups.txt"), "w") as output_file:
        for seed, members in ordered:
            total = sum(c for _, _, c in members)
            output_file.write(f"{seed} {total} {total / total_records:.2%}\n")
        if below_indexes:
            output_file.write(f"# {below}\n")


def _write_discovery(ordered: list, total_records: int, handle) -> None:
    """Stream one line per grouped index to ``handle``."""
    for seed, members in ordered:
        for dist, index, cnt in members:
            handle.write(f"{seed} {dist} {index} {cnt} {cnt / total_records:.2%}\n")


def _run_patterns(
    path: pathlib.Path,
    patterns_list: list[str],
    args: argparse.Namespace,
    out_dir: pathlib.Path | None,
) -> None:
    """Search one input file for each pattern."""
    for i, pattern in enumerate(patterns_list, start=1):
        logging.info(f"Processing pattern {i}: {pattern}")

        regex_pattern = compile_index_regex(pattern, args.distance, args.error_type)

        results = defaultdict(
            Counter, [(str(x), Counter()) for x in range(args.distance + 1)]
        )
        total_records = 0

        rust_kernel = None
        if args.engine == "rust":
            if biomate_kernels is None:
                logging.warning(
                    "biomate_kernels is not installed; falling back to the regex engine"
                )
            elif not re.fullmatch(r"[A-Za-z+]+", pattern):
                # The kernel matches letters literally: "ACG[TA]" or "ACG.A"
                # would give different results from the regex engine
                logging.warning(
                    f"Pattern '{pattern}' uses regex syntax; using the regex engine"
                )
            else:
                rust_kernel = biomate_kernels.IndexPattern(
                    pattern, args.distance, args.error_type
                )
                logging.info("Using the rust scan engine (biomate_kernels)")

        def count_match(seq: str) -> None:
            match = regex_pattern.search(seq)
            if match:
                fuzziness = str(sum(match.fuzzy_counts))
                results[fuzziness][(match.groups()[1:])] += 1

        def flush_chunk(chunk: list[str]) -> None:
            """Count one chunk with the rust engine; re-run fallbacks via regex."""
            aggregated, fallback = rust_kernel.scan_batch(chunk)
            for fuzziness, groups in aggregated.items():
                counter = results[fuzziness]
                for key, count in groups.items():
                    counter[tuple(key.split(" "))] += count
            for idx in fallback:
                count_match(chunk[idx])

        with dnaio.open(path) as reader:
            chunk: list[str] = []
            for record in reader:
                if rust_kernel is None:
                    count_match(record.sequence)
                else:
                    chunk.append(record.sequence)
                    if len(chunk) >= _KERNEL_CHUNK:
                        flush_chunk(chunk)
                        chunk = []
                total_records += 1
            if chunk:
                flush_chunk(chunk)
        logging.debug(f"Total records processed: {total_records}")
        if not [x for x in results.values() if x]:
            logging.warning("   No matches found.")
            continue

        if out_dir is not None:
            write_results(results, total_records, out_dir, file_index=i)
        else:
            for key, value in results.items():
                for seq, cnt in _sorted_by_count(value):
                    # print, not logging: without --output the log is quiet
                    print(
                        f"{key} {' '.join(list(seq))} {cnt} {cnt / total_records:.2%}"
                    )


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if not args.output:
        args.quiet = True
    setup_logging(args)

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
    if not patterns_list:
        logging.info("No index pattern given: grouping the header indexes")

    multiple = len(args.input) > 1
    jobs = []
    for path in args.input:
        out_dir = args.output
        if multiple and out_dir is not None:
            out_dir = out_dir / _output_name(path)
        jobs.append((path, patterns_list, args, out_dir))

    # Workers would interleave their stdout, so only files run in parallel
    if args.threads > 1 and multiple and args.output is not None:
        with Pool(min(args.threads, len(jobs))) as pool:
            pool.starmap(_process_file, jobs)
        return
    for job in jobs:
        if multiple and args.output is None:
            print(f"# {job[0]}")
        _process_file(*job)


def _process_file(
    path: pathlib.Path,
    patterns_list: list[str],
    args: argparse.Namespace,
    out_dir: pathlib.Path | None,
) -> None:
    """Search or discover the indexes of one input file."""
    logging.info(f"Input file: {path}")
    logging.info(f"Input type: {_detect_input_type(path.name)}")
    if patterns_list:
        _run_patterns(path, patterns_list, args, out_dir)
    else:
        _run_discovery(path, args, out_dir)
