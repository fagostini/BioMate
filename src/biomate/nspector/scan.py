"""FASTQ scanning for nspector: N positions, per-tile totals and file metadata."""

from __future__ import annotations

import dataclasses
import logging
import pathlib
import re
from collections import Counter

import dnaio

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
    """N-content scan of one FASTQ file.

    ``lanes``/``tiles``/``xs``/``ys``/``cycles`` are parallel lists over the
    N-containing reads; ``cycles`` holds 1-based N positions. ``tile_totals``
    maps (lane, tile) to the number of reads seen on that tile.
    """

    path: pathlib.Path
    sample: str
    read: str
    format_error: int | None = None
    lanes: list[int] = dataclasses.field(default_factory=list)
    tiles: list[int] = dataclasses.field(default_factory=list)
    xs: list[int] = dataclasses.field(default_factory=list)
    ys: list[int] = dataclasses.field(default_factory=list)
    cycles: list[list[int]] = dataclasses.field(default_factory=list)
    tile_totals: dict[tuple[int, int], int] = dataclasses.field(default_factory=dict)
    max_len: int = 0

    @property
    def stem(self) -> str:
        """File name without FASTQ suffixes."""
        return fastq_stem(self.path)

    @property
    def total_reads(self) -> int:
        """Number of reads scanned."""
        return sum(self.tile_totals.values())


def scan_file_python(path: pathlib.Path, threads: int = 0) -> FileScan:
    """Scan one FASTQ with dnaio (reference implementation)."""
    sample, read = parse_fastq_name(path)
    scan = FileScan(path=path, sample=sample, read=read)
    totals: Counter[tuple[int, int]] = Counter()

    with dnaio.open(path, open_threads=threads) as reader:
        for idx, record in enumerate(reader):
            coords = extract_tile_coordinates(record.name)
            if coords is None:
                scan.format_error = idx
                break
            lane, tile, x, y = coords
            totals[(lane, tile)] += 1
            sequence = record.sequence
            if len(sequence) > scan.max_len:
                scan.max_len = len(sequence)
            if "N" in sequence:
                scan.lanes.append(lane)
                scan.tiles.append(tile)
                scan.xs.append(x)
                scan.ys.append(y)
                scan.cycles.append([i + 1 for i, c in enumerate(sequence) if c == "N"])

    scan.tile_totals = dict(sorted(totals.items()))
    return scan


def _from_kernel(path: pathlib.Path, result: tuple) -> FileScan:
    """Convert one ``biomate_kernels.scan_n_content`` tuple to a FileScan."""
    fmt_err, lanes, tiles, xs, ys, cycles, totals, max_len = result
    sample, read = parse_fastq_name(path)
    return FileScan(
        path=path,
        sample=sample,
        read=read,
        format_error=fmt_err,
        lanes=list(lanes),
        tiles=list(tiles),
        xs=list(xs),
        ys=list(ys),
        cycles=[list(c) for c in cycles],
        tile_totals={(lane, tile): count for lane, tile, count in totals},
        max_len=max_len,
    )


def scan_files(
    paths: list[pathlib.Path], engine: str = "python", threads: int = 0
) -> list[FileScan]:
    """Scan every input file with the requested engine.

    The Rust engine falls back to Python when ``biomate_kernels`` is missing.
    """
    if engine == "rust":
        if biomate_kernels is None:
            logging.warning(
                "biomate_kernels is not available; falling back to the Python scan."
            )
        else:
            results = biomate_kernels.scan_n_content([str(p) for p in paths], threads)
            return [_from_kernel(p, r) for p, r in zip(paths, results)]
    return [scan_file_python(p, threads) for p in paths]
