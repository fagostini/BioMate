"""FASTQ scanning for nspector: N positions, per-tile totals and file metadata."""

from __future__ import annotations

import dataclasses
import logging
import os
import pathlib
import re
from collections import Counter, defaultdict
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor

import dnaio
import numpy
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
# Cycles at each end of the read used to classify where the first N falls.
EDGE_CYCLES = 10


def first_n_region(first: int, length: int) -> int:
    """Return 0 (first EDGE_CYCLES cycles), 2 (last EDGE_CYCLES) or 1 (middle).

    ``first`` is the 0-based index of the first N; start wins on short reads.
    """
    if first < EDGE_CYCLES:
        return 0
    if first + EDGE_CYCLES >= length:
        return 2
    return 1


def longest_run(cycles: list[int]) -> int:
    """Return the longest run of consecutive cycles in a sorted cycle list."""
    best = run = 1
    for prev, cur in zip(cycles, cycles[1:]):
        run = run + 1 if cur == prev + 1 else 1
        best = max(best, run)
    return best


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


UNDETERMINED = "Undetermined"


def project_of(path: pathlib.Path) -> str:
    """Return the project of a FASTQ from its folder layout.

    ``Undetermined`` files form their own project. Otherwise the project is the
    parent folder (BCL Convert ``<Sample_Project>/<file>``), or the folder
    above a ``Sample_*`` folder (bcl2fastq ``<Project>/Sample_<id>/<file>``).
    """
    if parse_fastq_name(path)[0] == UNDETERMINED:
        return UNDETERMINED
    parent = path.parent
    if parent.name.startswith("Sample_"):
        parent = parent.parent
    return parent.name or "."


def discover_fastqs(inputs: list[pathlib.Path]) -> list[pathlib.Path]:
    """Expand directories into the FASTQ files below them (recursive, sorted).

    Files are kept as given; duplicates (same resolved path) are dropped,
    keeping the first occurrence.
    """
    found: dict[pathlib.Path, pathlib.Path] = {}
    for item in inputs:
        candidates = (
            sorted(
                p
                for p in item.rglob("*")
                if p.is_file() and p.name.endswith(_FASTQ_SUFFIXES)
            )
            if item.is_dir()
            else [item]
        )
        for path in candidates:
            found.setdefault(path.resolve(), path)
    return list(found.values())


@dataclasses.dataclass
class FileScan:
    """Aggregated N-content scan of one FASTQ file.

    Counter lists are sorted and hold only non-zero entries; their size depends
    on lanes x tiles x cycles, never on the number of reads:

    * ``tiles``: (lane, tile, total, reads_with_n, n_bases, lowq_reads, bases),
      ``bases`` being the sum of read lengths
    * ``tile_cycle``: (lane, tile, cycle, reads_with_n, reads_with_lowq),
      1-based cycles, rows where either count is non-zero
    * ``n_per_read``: (lane, n_count, reads)
    * ``spatial``: (lane, tile, x // bin_size, y // bin_size,
      (cycle - 1) // cycle_window, n_bases)
    * ``n_shape``: (lane, longest_n_run, first_n_region, all_n, reads), see
      :func:`first_n_region`

    Low-quality ("lowq") counts only called bases (not N) whose Phred score
    is below ``min_qual``.
    """

    path: pathlib.Path
    sample: str
    read: str
    project: str = ""
    format_error: int | None = None
    max_len: int = 0
    tiles: list[tuple[int, int, int, int, int, int, int]] = dataclasses.field(
        default_factory=list
    )
    tile_cycle: list[tuple[int, int, int, int, int]] = dataclasses.field(
        default_factory=list
    )
    n_per_read: list[tuple[int, int, int]] = dataclasses.field(default_factory=list)
    spatial: list[tuple[int, int, int, int, int, int]] = dataclasses.field(
        default_factory=list
    )
    n_shape: list[tuple[int, int, int, int, int]] = dataclasses.field(
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


class LowQCounter:
    """Count called (non-N) bases below ``min_qual`` per tile and cycle.

    Reads are buffered per (lane, tile, length) and compared in numpy
    batches; a per-base Python loop made the Python engine ~13x slower.
    At most ``batch_size`` reads are buffered in total, whatever the number
    of tiles.
    """

    def __init__(self, min_qual: int, batch_size: int = 8192):
        """Set the Phred+33 threshold; ``min_qual`` 0 counts nothing."""
        self.threshold = min_qual + 33
        self.batch_size = batch_size
        self.buffers: dict[tuple[int, int, int], tuple[list, list]] = defaultdict(
            lambda: ([], [])
        )
        self.cycles: dict[tuple[int, int], numpy.ndarray] = {}
        self.reads: Counter[tuple[int, int]] = Counter()
        self.buffered = 0

    def add(self, lane: int, tile: int, sequence: str, qualities: str) -> None:
        """Buffer one read, counting the buffer when it is full."""
        if self.threshold <= 33 or not sequence:
            return
        key = (lane, tile, len(sequence))
        seqs, quals = self.buffers[key]
        seqs.append(sequence)
        quals.append(qualities)
        self.buffered += 1
        if self.buffered >= self.batch_size:
            self._flush_all()

    def _flush(self, key: tuple[int, int, int]) -> None:
        lane, tile, length = key
        seqs, quals = self.buffers.pop(key)
        seq = numpy.frombuffer("".join(seqs).encode("ascii"), numpy.uint8)
        qual = numpy.frombuffer("".join(quals).encode("ascii"), numpy.uint8)
        low = ((qual < self.threshold) & (seq != ord("N"))).reshape(-1, length)
        counts = self.cycles.get((lane, tile))
        if counts is None or len(counts) < length:
            grown = numpy.zeros(length, numpy.int64)
            if counts is not None:
                grown[: len(counts)] = counts
            counts = self.cycles[(lane, tile)] = grown
        counts[:length] += low.sum(axis=0)
        self.reads[(lane, tile)] += int(low.any(axis=1).sum())

    def _flush_all(self) -> None:
        for key in list(self.buffers):
            self._flush(key)
        self.buffered = 0

    def finish(self) -> Counter[tuple[int, int, int]]:
        """Count what is left and return {(lane, tile, cycle): reads with low Q}."""
        self._flush_all()
        return Counter(
            {
                (lane, tile, cycle + 1): int(c)
                for (lane, tile), counts in self.cycles.items()
                for cycle, c in enumerate(counts)
                if c
            }
        )


def scan_file_python(
    path: pathlib.Path,
    threads: int = 0,
    bin_size: int = 1000,
    raw: RawWriter | None = None,
    min_qual: int = 20,
    cycle_window: int = 10,
) -> FileScan:
    """Scan one FASTQ with dnaio (reference implementation of the Rust kernel)."""
    sample, read = parse_fastq_name(path)
    scan = FileScan(path=path, sample=sample, read=read, project=project_of(path))
    tiles: dict[tuple[int, int], list[int]] = {}
    tile_cycle: Counter[tuple[int, int, int]] = Counter()
    lowq = LowQCounter(min_qual)
    n_per_read: Counter[tuple[int, int]] = Counter()
    spatial: Counter[tuple[int, int, int, int, int]] = Counter()
    n_shape: Counter[tuple[int, int, int, int]] = Counter()

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
            acc = tiles.setdefault((lane, tile), [0, 0, 0, 0])
            acc[0] += 1
            acc[3] += len(sequence)
            lowq.add(lane, tile, sequence, record.qualities)
            if "N" not in sequence:
                continue
            cycles = [i + 1 for i, c in enumerate(sequence) if c == "N"]
            acc[1] += 1
            acc[2] += len(cycles)
            for cycle in cycles:
                tile_cycle[(lane, tile, cycle)] += 1
            n_per_read[(lane, len(cycles))] += 1
            all_n = int(len(cycles) == len(sequence))
            region = first_n_region(cycles[0] - 1, len(sequence))
            n_shape[(lane, longest_run(cycles), region, all_n)] += 1
            xb, yb = x // bin_size, y // bin_size
            for cycle in cycles:
                spatial[(lane, tile, xb, yb, (cycle - 1) // cycle_window)] += 1
            if raw is not None:
                raw.add(lane, tile, x, y, cycles)

    if raw is not None:
        raw.flush()
    lowq_cycle = lowq.finish()
    scan.tiles = sorted(
        (lane, tile, *acc[:3], lowq.reads[(lane, tile)], acc[3])
        for (lane, tile), acc in tiles.items()
    )
    scan.tile_cycle = sorted(
        (*k, tile_cycle.get(k, 0), lowq_cycle.get(k, 0))
        for k in tile_cycle.keys() | lowq_cycle.keys()
    )
    scan.n_per_read = sorted((*k, v) for k, v in n_per_read.items())
    scan.spatial = sorted((*k, v) for k, v in spatial.items())
    scan.n_shape = sorted((*k, v) for k, v in n_shape.items())
    return scan


def _from_kernel(path: pathlib.Path, result: tuple) -> FileScan:
    """Convert one ``biomate_kernels.scan_n_content`` tuple to a FileScan."""
    fmt_err, max_len, tiles, tile_cycle, n_per_read, spatial, n_shape = result
    sample, read = parse_fastq_name(path)
    return FileScan(
        path=path,
        sample=sample,
        read=read,
        project=project_of(path),
        format_error=fmt_err,
        max_len=max_len,
        tiles=[tuple(t) for t in tiles],
        tile_cycle=[tuple(t) for t in tile_cycle],
        n_per_read=[tuple(t) for t in n_per_read],
        spatial=[tuple(t) for t in spatial],
        n_shape=[tuple(t) for t in n_shape],
    )


def resolve_engine(engine: str) -> str:
    """Map ``auto`` to ``rust`` when biomate_kernels imports, else ``python``."""
    if engine == "auto":
        return "rust" if biomate_kernels is not None else "python"
    if engine == "rust" and biomate_kernels is None:
        raise RuntimeError(
            "--engine rust needs biomate_kernels, which is not installed in this "
            "environment. Build it (cd rust/biomate_kernels && maturin develop "
            "--release) or use --engine auto/python."
        )
    return engine


def scan_files(
    paths: list[pathlib.Path],
    engine: str = "auto",
    threads: int = 0,
    bin_size: int = 1000,
    raw_dir: pathlib.Path | None = None,
    min_qual: int = 20,
    cycle_window: int = 10,
    done: Callable[[pathlib.Path], None] | None = None,
) -> list[FileScan]:
    """Scan every input file with the requested engine.

    With ``raw_dir`` set, per-read records are streamed to
    ``raw_dir/<stem>_n_sequences/`` and the Python engine is used.
    ``done(path)`` is called as each file finishes (e.g. to advance a
    progress bar); with the Rust engine it runs on a worker thread.
    """
    engine = resolve_engine(engine)

    def finished(scan: FileScan) -> FileScan:
        if done is not None:
            done(scan.path)
        return scan

    if raw_dir is not None:
        # shortcut: raw export only in the Python engine; add a Rust writer if --raw is used on big runs
        if engine == "rust":
            logging.info("--raw uses the Python scan engine.")
        return [
            finished(
                scan_file_python(
                    p,
                    threads,
                    bin_size,
                    RawWriter(raw_dir / f"{fastq_stem(p)}_n_sequences"),
                    min_qual,
                    cycle_window,
                )
            )
            for p in paths
        ]
    if engine == "rust":
        # 0 = one file per available core (sched_getaffinity respects taskset/cgroups)
        workers = threads or len(os.sched_getaffinity(0))

        def scan_one(path: pathlib.Path) -> FileScan:
            # one file per call so progress can be reported; the kernel
            # releases the GIL, so the threads scan in parallel
            result = biomate_kernels.scan_n_content(
                [str(path)], 1, bin_size, min_qual, cycle_window
            )
            return finished(_from_kernel(path, result[0]))

        with ThreadPoolExecutor(max_workers=min(workers, max(len(paths), 1))) as pool:
            return list(pool.map(scan_one, paths))
    return [
        finished(
            scan_file_python(
                p, threads, bin_size, min_qual=min_qual, cycle_window=cycle_window
            )
        )
        for p in paths
    ]
