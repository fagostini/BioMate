"""Waif main script.

Report barcodes that remain undetermined in every demultiplexing run of a
lane.

For each lane that appears in two or more datasets (bcl-convert runs stored in
``<input_path>/<Dataset>/Stats/Stats.json`` directories), the barcodes that
are undetermined in *all* of those datasets are reported with an estimated
read count: the mean of the datasets' counts, rounded up to a multiple of 10.
"""

import argparse
import gzip
import hashlib
import io
import json
import logging
import os
import pickle
import queue
import sys
import threading
import zlib
from collections import Counter, defaultdict
from collections.abc import Iterable
from concurrent.futures import ProcessPoolExecutor

import pathlib

import polars

from biomate.setup import setup_logging

STATS_FILE_NAME = "Stats.json"
EXCLUDED_DIR = "Demultiplexing"
INDEX_SEP = "+"
COUNT_ROUND_TO = 10
COUNT_COLUMN_PREFIX = "Count_"
THREADED_MIN_SIZE = 8 << 20
_INFLATE_CHUNK = 16 << 20
_INDEX_READ_CHUNK = 1 << 20


def init_parser(subparsers: argparse._SubParsersAction) -> argparse.ArgumentParser:
    """Initialise module subparser."""
    parser = subparsers.add_parser(
        __name__.split(".")[-1],
        description=(
            "Report barcodes that remain undetermined in every demultiplexing "
            "run of a lane."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help=(
            "Report barcodes that remain undetermined in every demultiplexing "
            "run of a lane."
        ),
    )
    parser.add_argument(
        "--input-path",
        type=pathlib.Path,
        required=True,
        help="pathlib.Path to the input directory containing the stats files.",
    )
    parser.add_argument(
        "--output-file",
        type=pathlib.Path,
        required=False,
        help=(
            "File to write the results to. If not specified, the results are "
            "written to stdout as CSV."
        ),
    )
    parser.add_argument(
        "--exact",
        action="store_true",
        help=(
            "Enable exact counting of undetermined indexes by reading the "
            "undetermined FASTQ files."
        ),
    )
    parser.add_argument(
        "--sample-sheet",
        type=pathlib.Path,
        required=False,
        help=(
            "pathlib.Path to the Sample Sheet file (required for --exact). If not "
            "provided, the module looks for SampleSheet.csv in the input "
            "directory."
        ),
    )
    parser.add_argument(
        "--block-cache",
        type=pathlib.Path,
        default=None,
        metavar="DIR",
        help=(
            "Directory for a small gzip member-offset index of the "
            "undetermined files; enables forked per-file workers that "
            "decompress and count in parallel (effective when there are "
            "fewer files than cores). Each file is scanned once to locate "
            "member boundaries (a few minutes per large file); the index "
            "file is tiny and the source is never copied or re-compressed."
        ),
    )
    parser.set_defaults(parse=validate_args, run=main)

    return parser


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """Validate the command line arguments."""
    if not args.input_path.is_dir():
        raise argparse.ArgumentTypeError(
            f"Input path {args.input_path} does not exist or is not a directory."
        )
    if args.output_file is not None and not args.output_file.parent.is_dir():
        logging.info(
            f"Output directory {args.output_file.parent} does not exist. Creating it..."
        )
        args.output_file.parent.mkdir(parents=True, exist_ok=True)
    if args.exact:
        if args.sample_sheet is not None and not args.sample_sheet.is_file():
            raise argparse.ArgumentTypeError(
                f"Sample Sheet {args.sample_sheet} does not exist."
            )
        sample_sheet = args.sample_sheet or args.input_path / "SampleSheet.csv"
        if not sample_sheet.is_file():
            raise argparse.ArgumentTypeError(
                f"Sample Sheet not found in {args.input_path} "
                "(required for --exact); provide one with --sample-sheet."
            )
    if args.block_cache is not None:
        if sys.platform != "linux":
            raise argparse.ArgumentTypeError(
                "--block-cache is only supported on Linux (it relies on os.fork)."
            )
        args.block_cache.mkdir(parents=True, exist_ok=True)
    return args


def get_shared_lanes_from_samplesheet(
    samplesheet: pathlib.Path,
) -> list[int]:
    """Return the lanes carrying samples from two or more Sample_Projects.

    Rows without any index (``index``/``index2`` both empty) and rows with an
    empty project are ignored.
    """
    projects_by_lane: defaultdict[str, set[str]] = defaultdict(set)
    spreadsheet_version = None
    with open(samplesheet, "r") as f:
        for line in f:
            line = line.strip()
            if spreadsheet_version:
                if line:
                    parts = line.split(",")
                    if spreadsheet_version == 1:
                        lane, index1, index2, project = (
                            parts[1],
                            parts[5],
                            parts[6],
                            parts[11] if len(parts) > 11 else "",
                        )
                    else:
                        lane, index1, index2, project = (
                            parts[0],
                            parts[3],
                            parts[4],
                            parts[5] if len(parts) > 5 else "",
                        )
                    if (index1 or index2) and project:
                        projects_by_lane[lane].add(project)
            elif line.startswith("Lane,"):
                spreadsheet_version = 2
            elif line.startswith("FCID,"):
                spreadsheet_version = 1
            else:
                continue

    return sorted(
        int(lane) for lane, projects in projects_by_lane.items() if len(projects) > 1
    )


def _count_header_lines(lines: Iterable[str]) -> Counter[str]:
    """Count the index pairs of the header lines of a FASTQ stream."""
    counts: Counter[str] = Counter()
    for line in lines:
        if line.startswith("@"):
            counts[line.rsplit(":", 1)[-1].rstrip()] += 1
    return counts


def _count_serial(f: pathlib.Path) -> Counter[str]:
    """One fused streaming pass: decompress and count line by line."""
    with gzip.open(f, "rt") as fh:
        return _count_header_lines(fh)


class _QueueReader(io.RawIOBase):
    """Read-only raw file object fed with decompressed chunks from a queue.

    ``readinto`` consumes each queued chunk via a moving index (no repeated
    copying of the remainder) and keeps returning 0 after EOF, as
    ``TextIOWrapper`` re-reads after the final line. A ``BaseException``
    queued by the producer is recorded in ``error``.
    """

    def __init__(self, q: queue.Queue) -> None:
        """Feed this reader with the chunks queued by ``q``."""
        super().__init__()
        self.q = q
        self.buf = b""
        self.pos = 0
        self.eof = False
        self.error: BaseException | None = None

    def readable(self) -> bool:
        """True while the stream is open."""
        return not self.closed

    def readinto(self, buf: bytearray) -> int:
        """Fill ``buf`` from the queued chunks; return the bytes read."""
        if self.eof:
            return 0
        while self.pos >= len(self.buf):
            item = self.q.get()
            if item is None:
                self.eof = True
                return 0
            if isinstance(item, BaseException):
                self.error = item
                self.eof = True
                return 0
            self.buf = item
            self.pos = 0
        n = min(len(buf), len(self.buf) - self.pos)
        buf[:n] = self.buf[self.pos : self.pos + n]
        self.pos += n
        return n


def _inflate_stream(f: pathlib.Path, q: queue.Queue, chunk: int | None = None) -> None:
    """Producer thread: stream decompressed data of a gzip file into ``q``.

    The per-chunk ``zlib`` inflate C call releases the GIL, so this thread
    runs truly in parallel with the counting consumer. Multi-member files
    are handled by re-creating the decompressor: on ``unused_data`` directly,
    or by reading the next chunk when a member ends exactly at a chunk edge
    (only a read returning EOF ends the stream, so no member is skipped).
    On failure the exception itself is queued (re-raised by the caller); a
    ``None`` sentinel terminates the stream in all cases.
    """
    size = _INFLATE_CHUNK if chunk is None else chunk
    try:
        d = zlib.decompressobj(47)
        with open(f, "rb") as fh:
            comp = fh.read(size)
            while comp:
                out = d.decompress(comp)
                if out:
                    q.put(out)
                if d.eof:
                    tail = d.flush()
                    if tail:
                        q.put(tail)
                    if d.unused_data:
                        comp = d.unused_data
                        d = zlib.decompressobj(47)
                        continue
                    comp = fh.read(size)
                    if not comp:
                        break
                    d = zlib.decompressobj(47)
                    continue
                comp = fh.read(size)
    except BaseException as e:  # noqa: BLE001 - re-raised by the caller
        q.put(e)
    finally:
        q.put(None)


def _count_threaded(f: pathlib.Path) -> Counter[str]:
    """Count with zlib inflate (producer thread) overlapped with counting."""
    q: queue.Queue = queue.Queue()
    reader = _QueueReader(q)
    th = threading.Thread(target=_inflate_stream, args=(f, q), daemon=True)
    th.start()
    counts = _count_header_lines(io.TextIOWrapper(reader, encoding="latin-1"))
    th.join()
    if reader.error is not None:
        raise reader.error
    return counts


def _count_undetermined(f: pathlib.Path, force_threaded: bool = False) -> Counter[str]:
    """Count the index pairs found in the headers of a gzipped FASTQ file.

    Small files use one fused streaming pass. For larger files the zlib
    inflate C call (which releases the GIL) runs in a producer thread,
    overlapped with the counting pass: ~15% faster than the fused pass in
    benchmarks (91MB file, 675MB decompressed: 2.04s -> 1.73s).
    """
    if force_threaded or f.stat().st_size >= THREADED_MIN_SIZE:
        return _count_threaded(f)
    return _count_serial(f)


def _member_index_for(f: pathlib.Path, cache_dir: pathlib.Path) -> pathlib.Path:
    """Return the cache path of the member index for ``f``."""
    digest = hashlib.sha256(str(f.resolve()).encode()).hexdigest()[:16]
    return cache_dir / f"{digest}_{f.name}.idx.json"


def _load_member_index(f: pathlib.Path, index_path: pathlib.Path) -> dict | None:
    """Return the member index if it still matches the source file."""
    try:
        index = json.loads(index_path.read_text())
    except (OSError, ValueError):
        return None
    st = f.stat()
    if (
        index.get("src") == str(f.resolve())
        and index.get("mtime_ns") == st.st_mtime_ns
        and index.get("size") == st.st_size
    ):
        return index
    return None


def _find_member_offsets(f: pathlib.Path) -> tuple[list[int], bool]:
    """Start byte offset of every gzip member, and whether they are line-safe.

    ``offsets[k]`` is where member ``k`` begins and ``offsets[-1]`` is the end
    of the file, so member ``k`` spans ``offsets[k]:offsets[k + 1]``. The
    boolean is True only when every member boundary falls on a FASTQ line
    boundary (each non-final member's decompressed content is empty or ends on
    a newline); only then can a byte range be split at member boundaries and
    counted line by line. One sequential decompress pass; nothing is retained.
    A misaligned boundary aborts early, so a bad file is rejected in seconds.
    """
    offsets = [0]
    size = f.stat().st_size
    d = zlib.decompressobj(47)
    pending = b""
    comp_read = 0
    last_byte: int | None = None  # last decompressed byte of the current member
    with open(f, "rb") as fh:
        while True:
            if not pending:
                block = fh.read(_INDEX_READ_CHUNK)
                comp_read += len(block)
                if not block:
                    break
                pending = block
            out = d.decompress(pending)
            if out:
                last_byte = out[-1]
            leftover = d.unused_data
            if d.eof:
                tail = d.flush()
                if tail:
                    last_byte = tail[-1]
                end = comp_read - len(leftover)
                offsets.append(end)
                if end < size and last_byte not in (None, 10):
                    return offsets, False
                d = zlib.decompressobj(47)
                last_byte = None
            pending = leftover
    return offsets, True


def _build_member_index(f: pathlib.Path, index_path: pathlib.Path) -> dict:
    """Record every gzip member's byte offset in ``index_path``.

    One sequential decompress pass - no re-compression, no large cache file.
    The index is validated against the source by path, mtime and size.
    """
    offsets, aligned = _find_member_offsets(f)
    index_path.parent.mkdir(parents=True, exist_ok=True)
    st = f.stat()
    index = {
        "src": str(f.resolve()),
        "mtime_ns": st.st_mtime_ns,
        "size": st.st_size,
        "aligned": aligned,
        "offsets": offsets,
    }
    tmp = index_path.with_name(f"{index_path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(index))
    os.replace(tmp, index_path)
    return index


class BoundedReader(io.RawIOBase):
    """Read-only view of ``path[start:end]`` for ``gzip.GzipFile``.

    Lets a worker decompress a contiguous run of gzip members (a byte range
    whose bounds are member boundaries) without loading the whole file.
    """

    def __init__(self, path: pathlib.Path, start: int, end: int) -> None:
        """Expose ``path[start:end]`` as a seekable byte stream."""
        super().__init__()
        self._fh = open(path, "rb")
        self._fh.seek(start)
        self._start = start
        self._end = end
        self._pos = start

    def readable(self) -> bool:
        """True while the underlying file is open."""
        return not self.closed

    def readinto(self, buf: bytearray) -> int:
        """Read up to ``len(buf)`` bytes from the bounded range."""
        n = min(len(buf), self._end - self._pos)
        if n <= 0:
            return 0
        data = self._fh.read(n)
        self._pos += len(data)
        buf[: len(data)] = data
        return len(data)

    def seekable(self) -> bool:
        """The range view is seekable."""
        return True

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        """Seek within the bounded range, clamped to its ends."""
        if whence == io.SEEK_SET:
            target = self._start + offset
        elif whence == io.SEEK_END:
            target = self._end + offset
        else:
            target = self._pos + offset
        target = max(self._start, min(target, self._end))
        self._pos = target
        return self._pos

    def tell(self) -> int:
        """Return the current position within the bounded range."""
        return self._pos

    def close(self) -> None:
        """Close the underlying file."""
        if self._fh is not None:
            self._fh.close()
            self._fh = None
        super().close()


def _count_bounded_range(path: pathlib.Path, start: int, end: int) -> Counter[bytes]:
    """Count header barcodes (bytes keys) in ``path[start:end]``.

    ``start``/``end`` are gzip member boundaries, so the range is a whole set
    of members; ``GzipFile`` decompresses them without loading the file.
    """
    counts: Counter[bytes] = Counter()
    with gzip.GzipFile(fileobj=BoundedReader(path, start, end), mode="rb") as fh:
        for line in fh:
            if line.startswith(b"@"):
                counts[line.rsplit(b":", 1)[-1].rstrip()] += 1
    return counts


def _index_child(
    path: pathlib.Path,
    start: int,
    end: int,
    pipes: list[tuple[int, int]],
    idx: int,
) -> None:
    """Forked worker: count ``path[start:end]``, send the Counter on a pipe."""
    my_w = pipes[idx][1]
    for r, w in pipes:
        os.close(r)
    for i, (r, w) in enumerate(pipes):
        if i != idx:
            os.close(w)
    try:
        payload = pickle.dumps(_count_bounded_range(path, start, end))
    except BaseException as e:  # noqa: BLE001 - sent to the parent
        payload = pickle.dumps(e)
    try:
        while payload:
            written = os.write(my_w, payload)
            payload = payload[written:]
    finally:
        os.close(my_w)
        # _exit: skip atexit/GC so inherited state is untouched.
        os._exit(0)


def _read_pipe_result(r: int) -> object:
    """Read and unpickle the payload written by a forked index worker."""
    chunks: list[bytes] = []
    while chunk := os.read(r, 1 << 20):
        chunks.append(chunk)
    blob = b"".join(chunks)
    if not blob:
        raise RuntimeError("index worker produced no result (killed?)")
    return pickle.loads(blob)


def _count_from_index(f: pathlib.Path, index: dict, budget: int) -> Counter[str]:
    """Count ``f`` with one forked worker per contiguous run of members."""
    offsets = index["offsets"]
    n_members = len(offsets) - 1
    n_workers = min(budget, n_members)
    if n_workers <= 1:
        counts: Counter[bytes] = _count_bounded_range(f, offsets[0], offsets[-1])
    else:
        pipes = [os.pipe() for _ in range(n_workers)]
        pids: list[int] = []
        results: list[object] = []
        try:
            for i in range(n_workers):
                lo = (i * n_members) // n_workers
                hi = ((i + 1) * n_members) // n_workers
                start, end = offsets[lo], offsets[hi]
                pid = os.fork()
                if pid == 0:
                    _index_child(f, start, end, pipes, i)
                    os._exit(1)  # unreachable: _index_child never returns
                pids.append(pid)
            for r, w in pipes:
                os.close(w)
            # Drain every pipe BEFORE waitpid: a child blocks on a full pipe
            # (64KB) until read, so waiting first would deadlock.
            for i in range(n_workers):
                results.append(_read_pipe_result(pipes[i][0]))
            for i, pid in enumerate(pids):
                _, status = os.waitpid(pid, 0)
                if not os.WIFEXITED(status) or os.WEXITSTATUS(status) != 0:
                    raise RuntimeError(
                        f"index worker {i} exited abnormally (status={status})"
                    )
        finally:
            for r, w in pipes:
                for fd in (r, w):
                    try:
                        os.close(fd)
                    except OSError:
                        pass
        counts = Counter()
        for r in results:
            if isinstance(r, BaseException):
                raise RuntimeError(f"index counting worker failed: {r!r}") from r
            counts.update(r)  # type: ignore[arg-type]
    # Workers count with bytes keys; decode once for the (str-keyed) result.
    return Counter({k.decode("latin-1"): v for k, v in counts.items()})


def _count_for_file(
    f: pathlib.Path, budget: int, cache_dir: pathlib.Path | None
) -> Counter[str]:
    """Pick the counting strategy for a file given its core budget."""
    if budget < 2 or f.stat().st_size < THREADED_MIN_SIZE:
        return _count_serial(f)
    if cache_dir is not None:
        index_path = _member_index_for(f, cache_dir)
        index = _load_member_index(f, index_path)
        if index is None:
            index = _build_member_index(f, index_path)
        # The parallel count is only exact when every member boundary is a
        # line boundary; otherwise a header could straddle two workers.
        if index.get("aligned") and len(index["offsets"]) - 1 >= 2:
            return _count_from_index(f, index, budget)
        # Single-member or line-misaligned: count serially (always correct).
        return _count_serial(f)
    return _count_threaded(f)


def _count_file(
    task: tuple[int, pathlib.Path, int, pathlib.Path | None],
) -> tuple[int, pathlib.Path, dict[str, int]]:
    """Process-pool worker: count one FASTQ file (module-level so it pickles)."""
    lane, f, budget, cache_dir = task
    return lane, f, dict(_count_for_file(f, budget, cache_dir))


def _barcode_frame(counts: Counter, lane: int, f: pathlib.Path) -> polars.DataFrame:
    """Build the long barcode table of one counted file."""
    return (
        polars.DataFrame(
            {"Barcode": list(counts), "Count": list(counts.values())},
            schema={"Barcode": polars.String, "Count": polars.Int64},
        )
        .with_columns(
            polars.lit(lane, dtype=polars.Int64).alias("Lane"),
            polars.lit(f.parent.name).alias("Dataset"),
            polars.col("Barcode").str.len_chars().alias("Barcode_len"),
        )
        .select(["Lane", "Barcode", "Count", "Dataset", "Barcode_len"])
        .with_columns(
            polars.col("Barcode")
            .str.split_exact(INDEX_SEP, 1)
            .struct.rename_fields(["index1", "index2"])
            .alias("indexes")
        )
        .unnest("indexes")
        # A single-index barcode (no separator) yields a null index2; normalize
        # it to "" so a dataset is either uniformly single (index2 == "") or
        # uniformly dual (index2 non-empty).
        .with_columns(polars.col("index2").fill_null(""))
    )


def fastqs_to_stats_frames(
    input_dir: pathlib.Path,
    selected_lanes: list[int],
    block_cache: pathlib.Path | None = None,
) -> polars.DataFrame:
    """Convert the selected undetermined FASTQ files to a long barcode table.

    Expects the FASTQ files directly under ``<Dataset>/`` (the dataset name is
    ``file.parent.name``); only R1 files are used, as the index information is
    in the R1 headers. Files below any ``Demultiplexing/`` directory are
    excluded (that is the primary bcl-convert demux output, not a re-run). If
    bcl-convert split a dataset into several 1M-read parts (``_001``, ``_002``,
    ...), their counts are summed, so each (Lane, Dataset, Barcode) row is
    unique. Per-file counting runs in a process pool (one worker per file, up
    to the core count); the result order matches the file order, so the
    output is deterministic. Each file gets a core budget of
    ``cpu_count // len(files)``: with enough slack and a ``block_cache``
    directory it is counted by forked workers over a gzip member-offset index
    (a single-member file falls back to the serial pass); with slack but no
    ``block_cache`` it uses the threaded inflate+count pipeline; otherwise a
    single fused serial pass.

    Raises:
        ValueError: If no undetermined FASTQ files are found for the
            selected lanes.
    """
    lanes = set(selected_lanes)
    files = [
        (lane, f)
        for lane in sorted(lanes)
        for f in sorted(input_dir.glob(f"**/Undetermined_*_L00{lane}_R1_*.fastq.gz"))
        if EXCLUDED_DIR not in f.parts
    ]
    if not files:
        raise ValueError(
            f"No undetermined FASTQ files found for lanes {sorted(lanes)} "
            f"under {input_dir}."
        )
    for lane in sorted(lanes - {lane for lane, _ in files}):
        logging.warning(f"No undetermined FASTQ files for lane {lane}; skipping")

    cpu = os.cpu_count() or 1
    budget = max(1, cpu // len(files))
    tasks = [(lane, f, budget, block_cache) for lane, f in files]
    workers = min(len(files), cpu)
    if workers > 1:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(_count_file, tasks, chunksize=1))
    else:
        results = [_count_file(task) for task in tasks]
    frames = [_barcode_frame(counts, lane, f) for lane, f, counts in results]
    return (
        polars.concat(frames)
        .group_by(["Lane", "Dataset", "Barcode"])
        .agg(
            polars.sum("Count").alias("Count"),
            polars.col("Barcode_len").first().alias("Barcode_len"),
            polars.col("index1").first().alias("index1"),
            polars.col("index2").first().alias("index2"),
        )
    )


def load_stats_frames(input_dir: pathlib.Path) -> polars.DataFrame:
    """Load every Stats.json under ``input_dir`` into a long barcode table.

    Expected layout: ``input_dir/<Dataset>/Stats/Stats.json``; the dataset
    name is the directory that contains the ``Stats`` folder
    (``file.parents[1]``). Files below any ``Demultiplexing/`` directory are
    excluded (that is the primary bcl-convert demux output, not a re-run).

    Raises:
        ValueError: If no Stats.json files are found, or if none of them
            carries an ``UnknownBarcodes`` field.
    """
    stats_files = [
        f
        for f in input_dir.glob(f"**/{STATS_FILE_NAME}")
        if EXCLUDED_DIR not in f.parts
    ]
    if not stats_files:
        raise ValueError(f"No {STATS_FILE_NAME} files found under {input_dir}.")

    frames = []
    for f in stats_files:
        try:
            frame = polars.read_json(f).select("UnknownBarcodes")
        except polars.exceptions.ColumnNotFoundError:
            logging.warning(f"{f} has no 'UnknownBarcodes' field; skipping")
            continue
        frame = (
            frame.explode("UnknownBarcodes", empty_as_null=True)
            .unnest("UnknownBarcodes")
            .unnest("Barcodes")
            .unpivot(index="Lane", variable_name="Barcode", value_name="Count")
            .drop_nulls()
            .with_columns(
                polars.lit(f.parents[1].name).alias("Dataset"),
                polars.col("Barcode").str.len_chars().alias("Barcode_len"),
            )
            .with_columns(
                polars.col("Barcode")
                .str.split_exact(INDEX_SEP, 1)
                .struct.rename_fields(["index1", "index2"])
                .alias("indexes")
            )
            .unnest("indexes")
            # A single-index barcode (no separator) yields a null index2;
            # normalize it to "" (see _barcode_frame).
            .with_columns(polars.col("index2").fill_null(""))
        )
        frames.append(frame)

    if not frames:
        raise ValueError(
            f"No usable {STATS_FILE_NAME} files found (all lacked 'UnknownBarcodes')."
        )
    return polars.concat(frames)


def validate_table(table: polars.DataFrame) -> None:
    """Fail fast on input that violates the expected data contract.

    Every barcode must be ``index1`` or ``index1+index2`` with a non-empty
    index1, and each (Dataset, Lane) must use a single index1 length and a
    single index2 length.

    Raises:
        ValueError: If the table violates the data contract.
    """
    # index2 may be empty (a single-index barcode); index1 must not be.
    bad = table.filter(polars.col("index1").is_null() | (polars.col("index1") == ""))
    if bad.height:
        examples = bad.select("Dataset", "Barcode").unique().head(5)
        raise ValueError(
            "Found barcodes with an empty index1; every barcode must be "
            f"'index1' or 'index1{INDEX_SEP}index2'. Examples:\n{examples}"
        )
    mixed = (
        table.group_by(["Dataset", "Lane"])
        .agg(
            polars.col("index1").str.len_chars().n_unique().alias("i1_lengths"),
            polars.col("index2").str.len_chars().n_unique().alias("i2_lengths"),
        )
        .filter((polars.col("i1_lengths") > 1) | (polars.col("i2_lengths") > 1))
    )
    if mixed.height:
        raise ValueError(
            f"Expected one i1/i2 index length per (Dataset, Lane), found:\n{mixed}"
        )


def find_shared_lanes(table: polars.DataFrame) -> list[int]:
    """Return the lanes that appear in more than one dataset."""
    return (
        table.unique(subset=["Dataset", "Lane"])
        .group_by("Lane")
        .agg(polars.len().alias("Count"))
        .filter(polars.col("Count") > 1)
        .get_column("Lane")
        .to_list()
    )


def index_lengths(frame: polars.DataFrame) -> tuple[int, int]:
    """Return (i1 length, i2 length) of ``frame`` (unique by data contract)."""
    return (
        len(frame.head(1).get_column("index1").to_list()[0]),
        len(frame.head(1).get_column("index2").to_list()[0]),
    )


def has_index2(frame: polars.DataFrame) -> bool:
    """True if ``frame`` carries a non-empty second index (a dual-index dataset).

    A dataset is uniformly single-index (index2 == "") or uniformly dual-index
    (index2 non-empty) by the ``validate_table`` data contract, so one row
    decides.
    """
    return index_lengths(frame)[1] > 0


def _join_single_index(
    running: polars.DataFrame,
    other: polars.DataFrame,
    ds: str,
    count_cols: list[str],
    running_i1: int,
    other_i1: int,
) -> tuple[polars.DataFrame, bool]:
    """Join the running intersection with a step where at least one side is
    single-index, comparing on index1 only.

    Trims the longer index1 to the shorter (either side), collapses both sides
    over (index1, Lane) — summing each dataset's counts and discarding index2 —
    and inner-joins on index1 alone. Returns the new running frame (index2 set
    to "" and Barcode set to index1) and False (the result carries no index2).
    """
    count_col = f"{COUNT_COLUMN_PREFIX}{ds}"
    target_i1 = min(running_i1, other_i1)
    if running_i1 > target_i1:
        running = running.with_columns(polars.col("index1").str.slice(0, target_i1))
    if other_i1 > target_i1:
        other = other.with_columns(polars.col("index1").str.slice(0, target_i1))
    running = (
        running.group_by(["index1", "Lane"])
        .agg(
            [polars.sum(c).alias(c) for c in count_cols]
            + [polars.col("Dataset").first().alias("Dataset")]
        )
        .select(["Lane", *count_cols, "Dataset", "index1"])
    )
    other = (
        other.group_by(["index1", "Lane"])
        .agg(polars.sum(count_col).alias(count_col))
        .select(["Lane", count_col, "index1"])
    )
    return (
        running.join(other, on=["index1", "Lane"], how="inner")
        .with_columns(
            polars.col("index1").alias("Barcode"),
            polars.lit("", dtype=polars.String).alias("index2"),
        )
        .select(
            [
                "Lane",
                "Barcode",
                *count_cols,
                count_col,
                "Dataset",
                "index1",
                "index2",
            ]
        ),
        False,
    )


def intersect_lane_datasets(lane_dt: polars.DataFrame) -> polars.DataFrame:
    """Compute the barcodes undetermined in *all* of a lane's datasets.

    ``lane_dt`` must be sorted so longer barcodes come first (see ``main``);
    with ``unique(maintain_order=True)`` the first dataset is therefore the
    longest-index one (the hub), so the side that gets truncated is
    deterministic. The running intersection is inner-joined against each
    subsequent dataset, truncating and summing counts when the next dataset
    has shorter indexes.

    Datasets may be dual-index (``index1+index2``) or single-index
    (``index1`` only, i.e. index2 == ""). When every dataset in the lane is
    dual-index, the join is on (index1, index2) and the result keeps both
    indexes. As soon as a single-index dataset joins, the comparison drops to
    index1 only: the longer index1 is trimmed to the shorter (whichever side),
    index2 is discarded, and the result is reported in the single-index form.

    Each dataset's raw count is carried in its own ``Count_<dataset>`` column;
    the reported ``Count`` is their mean, rounded up to a multiple of
    ``COUNT_ROUND_TO`` with a single rounding at the end. The reported
    ``Barcode`` comes from the last-joined (shortest) dataset, i.e. it is in
    the shortest index-length space of the intersection.

    Known limitation: two datasets with equal total index length but crossed
    per-index lengths (e.g. 10+8 vs 8+10) tie in the sort and the name order
    then decides the hub, which can leave one index untruncated and make that
    lane's intersection come out empty.
    """
    lane_datasets = lane_dt["Dataset"].unique(maintain_order=True).to_list()
    by_dataset = {
        ds: lane_dt.filter(polars.col("Dataset") == ds) for ds in lane_datasets
    }

    hub = lane_datasets[0]
    running = by_dataset[hub].rename({"Count": f"{COUNT_COLUMN_PREFIX}{hub}"})
    running_has_i2 = has_index2(running)

    for ds in lane_datasets[1:]:
        other = by_dataset[ds].rename({"Count": f"{COUNT_COLUMN_PREFIX}{ds}"})
        count_cols = [c for c in running.columns if c.startswith(COUNT_COLUMN_PREFIX)]
        running_i1, running_i2 = index_lengths(running)
        other_i1, other_i2 = index_lengths(other)

        if not (running_has_i2 and has_index2(other)):
            # A single-index dataset is present (or the running intersection is
            # already single-index): compare on index1 only and drop index2.
            running, running_has_i2 = _join_single_index(
                running, other, ds, count_cols, running_i1, other_i1
            )
            if running.is_empty():
                break
            continue

        # Both dual-index: the original all-dual path, unchanged.
        if running_i1 > other_i1 or running_i2 > other_i2:
            if running_i1 > other_i1:
                running = running.with_columns(
                    polars.col("index1").str.slice(0, other_i1)
                )
            if running_i2 > other_i2:
                running = running.with_columns(
                    polars.col("index2").str.slice(0, other_i2)
                )
            # Truncation can collapse distinct full barcodes to the same
            # (index1, index2); sum each dataset's counts over the collapse.
            running = (
                running.group_by(["index1", "index2", "Lane"])
                .agg(
                    [polars.sum(c).alias(c) for c in count_cols]
                    + [
                        polars.col("Barcode").first().alias("Barcode"),
                        polars.col("Dataset").first().alias("Dataset"),
                    ]
                )
                .select(["Lane", "Barcode", *count_cols, "Dataset", "index1", "index2"])
            )

        running = (
            running.join(other, on=["index1", "index2", "Lane"], how="inner")
            .select(
                [
                    "Lane",
                    "Barcode_right",
                    *count_cols,
                    f"{COUNT_COLUMN_PREFIX}{ds}",
                    "Dataset",
                    "index1",
                    "index2",
                ]
            )
            .rename({"Barcode_right": "Barcode"})
        )
        if running.is_empty():
            break

    count_cols = [c for c in running.columns if c.startswith(COUNT_COLUMN_PREFIX)]
    mean_rounded_up = (
        polars.sum_horizontal(count_cols) / (COUNT_ROUND_TO * len(count_cols))
    ).ceil() * COUNT_ROUND_TO
    running = running.with_columns(mean_rounded_up.cast(polars.Int64).alias("Count"))
    return running.select(["Lane", "Barcode", "Count"])


def finalize(pair_frames: list[polars.DataFrame]) -> polars.DataFrame:
    """Concat the per-lane results and add the per-lane percentage column."""
    if not pair_frames:
        return polars.DataFrame(
            schema={
                "Lane": polars.Int64,
                "Barcode": polars.String,
                "Count": polars.Int64,
                "Total %": polars.Float64,
            }
        )
    return (
        polars.concat(pair_frames)
        .with_columns(
            (polars.col("Count") / polars.col("Count").sum() * 100)
            .over("Lane")
            .alias("Total %")
        )
        .sort(["Lane", "Count", "Barcode"], descending=[False, True, False])
    )


def main(args: argparse.Namespace) -> None:
    """Main function."""
    if args.output_file is None:
        args.quiet = True
    setup_logging(args)
    logging.info("Running waif module...")

    input_dir = args.input_path

    if args.exact:
        sample_sheet = args.sample_sheet or input_dir / "SampleSheet.csv"
        lanes_to_filter = get_shared_lanes_from_samplesheet(sample_sheet)
        if not lanes_to_filter:
            raise ValueError(
                "No lane with samples from two or more projects found in the "
                "Sample Sheet; nothing to do in --exact mode."
            )
        table = fastqs_to_stats_frames(input_dir, lanes_to_filter, args.block_cache)
    else:
        table = load_stats_frames(input_dir)
        lanes_to_filter = find_shared_lanes(table)

    present_lanes = set(table.get_column("Lane").unique().to_list())
    lanes_to_filter = [lane for lane in lanes_to_filter if lane in present_lanes]
    if not lanes_to_filter:
        message = (
            "No undetermined reads found for the selected lanes"
            if args.exact
            else "No lane appears in more than one dataset"
        )
        logging.info(f"{message}; nothing to report.")

    table = (
        table.filter(polars.col("Lane").is_in(lanes_to_filter))
        # Longer barcodes first: combined with unique(maintain_order=True) this
        # makes the longest-index dataset the hub (the truncated side).
        .sort(["Lane", "Barcode_len", "Dataset"], descending=[False, True, False])
        .drop(["Barcode_len"])
    )

    validate_table(table)

    pair_frames = [
        intersect_lane_datasets(table.filter(polars.col("Lane") == lane))
        for lane in lanes_to_filter
    ]

    results = finalize(pair_frames)

    if args.output_file is not None:
        results.write_csv(args.output_file)
        logging.info(f"Results saved to {args.output_file}.")
    else:
        sys.stdout.write(results.write_csv())
