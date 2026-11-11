"""FASTQ scanning for nspector: N positions, per-tile totals and file metadata."""

from __future__ import annotations

import dataclasses
import logging
import os
import pathlib
import re
from collections import Counter

import dnaio
import polars

try:
    import biomate_kernels
except ImportError:  # optional native kernel
    biomate_kernels = None


# Illumina BCL Convert / bcl2fastq file naming: Sample_S1_L001_R1_001.fastq.gz
_ILLUMINA_FASTQ_NAME = re.compile(
    r"^(?P<sample>.+?)_S\d+(?:_L(?P<lane>\d{3}))?_(?P<read>[RI]\d)_\d{3}$"
)
_FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")


def extract_tile_coordinates(read_name: str) -> tuple[int, int, int, int] | None:
    """Extract lane, tile and x, y coordinates from an Illumina read name.

    Illumina read name format: instrument:run:flowcell:lane:tile:x:y
    Returns (lane, tile, x, y) or None if the name format is invalid.
    """
    try:
        # Get the first part before any spaces (some formats have additional info)
        base_name = read_name.split(" ")[0]
        parts = base_name.split(":")
        if len(parts) < 7:
            return None
        return (int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6]))
    except (IndexError, ValueError):
        return None


def fastq_stem(path: pathlib.Path) -> str:
    """Return the file name without its FASTQ (and compression) suffix."""
    name = path.name
    for suffix in _FASTQ_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def parse_fastq_name(path: pathlib.Path) -> tuple[str, str]:
    """Return (sample, read) from an Illumina-style FASTQ file name.

    Falls back to (file stem, "R1") when the name does not follow the
    ``<sample>_S<n>[_L<lane>]_<R|I><n>_001`` convention.
    """
    stem = fastq_stem(path)
    match = _ILLUMINA_FASTQ_NAME.match(stem)
    if match is None:
        return stem, "R1"
    return match.group("sample"), match.group("read")


@dataclasses.dataclass
class FileScan:
    """Aggregated N-content scan of one FASTQ file.

    Counter lists are sorted and hold only non-zero entries; their size depends
    on lanes x tiles x cycles, never on the number of reads:

    * ``tiles``: (lane, tile, total, reads_with_n, n_bases)
    * ``tile_cycle``: (lane, tile, cycle, reads_with_n), 1-based cycles
    * ``n_per_read``: (lane, n_count, reads)
    * ``spatial``: (lane, tile, x // bin_size, y // bin_size, reads_with_n)
    """

    path: pathlib.Path
    sample: str
    read: str
    format_error: int | None = None
    max_len: int = 0
    tiles: list[tuple[int, int, int, int, int]] = dataclasses.field(
        default_factory=list
    )
    tile_cycle: list[tuple[int, int, int, int]] = dataclasses.field(
        default_factory=list
    )
    n_per_read: list[tuple[int, int, int]] = dataclasses.field(default_factory=list)
    spatial: list[tuple[int, int, int, int, int]] = dataclasses.field(
        default_factory=list
    )

    @property
    def stem(self) -> str:
        """File name without FASTQ suffixes."""
        return fastq_stem(self.path)

    @property
    def total_reads(self) -> int:
        """Number of reads scanned."""
        return sum(t[2] for t in self.tiles)

    @property
    def reads_with_n(self) -> int:
        """Number of reads with at least one N."""
        return sum(t[3] for t in self.tiles)


class RawWriter:
    """Stream per-read N records to ``<directory>/part-NNNNN.parquet`` in batches."""

    def __init__(self, directory: pathlib.Path, batch_size: int = 1_000_000):
        """Create the output directory and an empty batch."""
        self.directory = directory
        self.directory.mkdir(parents=True, exist_ok=True)
        self.batch_size = batch_size
        self.part = 0
        self.rows: dict[str, list] = {
            k: [] for k in ("lane", "tile", "x", "y", "cycles")
        }

    def add(self, lane: int, tile: int, x: int, y: int, cycles: list[int]) -> None:
        """Buffer one N-containing read, flushing full batches."""
        for key, value in zip(self.rows, (lane, tile, x, y, cycles)):
            self.rows[key].append(value)
        if len(self.rows["lane"]) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        """Write buffered rows as the next part file."""
        if not self.rows["lane"]:
            return
        polars.DataFrame(
            self.rows,
            schema={
                "lane": polars.UInt8,
                "tile": polars.UInt16,
                "x": polars.UInt32,
                "y": polars.UInt32,
                "cycles": polars.List(polars.UInt16),
            },
        ).write_parquet(self.directory / f"part-{self.part:05d}.parquet")
        self.part += 1
        self.rows = {k: [] for k in self.rows}


def scan_file_python(
    path: pathlib.Path,
    threads: int = 0,
    bin_size: int = 1000,
    raw: RawWriter | None = None,
) -> FileScan:
    """Scan one FASTQ with dnaio (reference implementation of the Rust kernel)."""
    sample, read = parse_fastq_name(path)
    scan = FileScan(path=path, sample=sample, read=read)
    tiles: dict[tuple[int, int], list[int]] = {}
    tile_cycle: Counter[tuple[int, int, int]] = Counter()
    n_per_read: Counter[tuple[int, int]] = Counter()
    spatial: Counter[tuple[int, int, int, int]] = Counter()

    with dnaio.open(path, open_threads=threads) as reader:
        for idx, record in enumerate(reader):
            coords = extract_tile_coordinates(record.name)
            if coords is None:
                scan.format_error = idx
                break
            lane, tile, x, y = coords
            sequence = record.sequence
            if len(sequence) > scan.max_len:
                scan.max_len = len(sequence)
            acc = tiles.setdefault((lane, tile), [0, 0, 0])
            acc[0] += 1
            if "N" not in sequence:
                continue
            cycles = [i + 1 for i, c in enumerate(sequence) if c == "N"]
            acc[1] += 1
            acc[2] += len(cycles)
            for cycle in cycles:
                tile_cycle[(lane, tile, cycle)] += 1
            n_per_read[(lane, len(cycles))] += 1
            spatial[(lane, tile, x // bin_size, y // bin_size)] += 1
            if raw is not None:
                raw.add(lane, tile, x, y, cycles)

    if raw is not None:
        raw.flush()
    scan.tiles = sorted((lane, tile, *acc) for (lane, tile), acc in tiles.items())
    scan.tile_cycle = sorted((*k, v) for k, v in tile_cycle.items())
    scan.n_per_read = sorted((*k, v) for k, v in n_per_read.items())
    scan.spatial = sorted((*k, v) for k, v in spatial.items())
    return scan


def _from_kernel(path: pathlib.Path, result: tuple) -> FileScan:
    """Convert one ``biomate_kernels.scan_n_content`` tuple to a FileScan."""
    fmt_err, max_len, tiles, tile_cycle, n_per_read, spatial = result
    sample, read = parse_fastq_name(path)
    return FileScan(
        path=path,
        sample=sample,
        read=read,
        format_error=fmt_err,
        max_len=max_len,
        tiles=[tuple(t) for t in tiles],
        tile_cycle=[tuple(t) for t in tile_cycle],
        n_per_read=[tuple(t) for t in n_per_read],
        spatial=[tuple(t) for t in spatial],
    )


def resolve_engine(engine: str) -> str:
    """Map ``auto`` to ``rust`` when biomate_kernels imports, else ``python``."""
    if engine == "auto":
        return "rust" if biomate_kernels is not None else "python"
    if engine == "rust" and biomate_kernels is None:
        logging.warning(
            "biomate_kernels is not available; falling back to the Python scan."
        )
        return "python"
    return engine


def scan_files(
    paths: list[pathlib.Path],
    engine: str = "auto",
    threads: int = 0,
    bin_size: int = 1000,
    raw_dir: pathlib.Path | None = None,
) -> list[FileScan]:
    """Scan every input file with the requested engine.

    With ``raw_dir`` set, per-read records are streamed to
    ``raw_dir/<stem>_n_sequences/`` and the Python engine is used.
    """
    engine = resolve_engine(engine)
    if raw_dir is not None:
        # shortcut: raw export only in the Python engine; add a Rust writer if --raw is used on big runs
        if engine == "rust":
            logging.info("--raw uses the Python scan engine.")
        return [
            scan_file_python(
                p,
                threads,
                bin_size,
                RawWriter(raw_dir / f"{fastq_stem(p)}_n_sequences"),
            )
            for p in paths
        ]
    if engine == "rust":
        # 0 = one file per available core (sched_getaffinity respects taskset/cgroups)
        workers = threads or len(os.sched_getaffinity(0))
        results = biomate_kernels.scan_n_content(
            [str(p) for p in paths], workers, bin_size
        )
        return [_from_kernel(p, r) for p, r in zip(paths, results)]
    return [scan_file_python(p, threads, bin_size) for p in paths]
