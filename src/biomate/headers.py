"""Count the indexes in FASTQ read headers (``... 1:N:0:ACGT+TTTT``)."""

import gzip
import hashlib
import io
import itertools
import json
import os
import pickle
import queue
import threading
import zlib
from collections import Counter
from collections.abc import Iterable

import pathlib

THREADED_MIN_SIZE = 8 << 20
_INFLATE_CHUNK = 16 << 20
_INDEX_READ_CHUNK = 1 << 20


def _count_header_lines(lines: Iterable[str]) -> Counter[str]:
    """Count the index pairs of the header lines of a FASTQ stream."""
    counts: Counter[str] = Counter()
    # Every 4th line is a header: a quality line may also start with "@" (Q31)
    for line in itertools.islice(lines, 0, None, 4):
        counts[line.rsplit(":", 1)[-1].rstrip()] += 1
    return counts


def _is_gzip(f: pathlib.Path) -> bool:
    """True when ``f`` starts with the gzip magic bytes."""
    with open(f, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _open_text(f: pathlib.Path) -> io.TextIOBase:
    """Open a gzipped or plain FASTQ file as text."""
    return gzip.open(f, "rt") if _is_gzip(f) else open(f)


def _count_serial(f: pathlib.Path) -> Counter[str]:
    """One fused streaming pass: decompress and count line by line."""
    with _open_text(f) as fh:
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
                # Cap every inflate call: a highly compressible chunk would
                # otherwise expand to gigabytes in one bytes object
                out = d.decompress(comp, _INFLATE_CHUNK)
                if out:
                    q.put(out)
                while not d.eof and (d.unconsumed_tail or len(out) == _INFLATE_CHUNK):
                    out = d.decompress(d.unconsumed_tail, _INFLATE_CHUNK)
                    if out:
                        q.put(out)
                if d.eof:
                    # Read before flush(): with a capped inflate the leftover is
                    # also in unconsumed_tail, and flush() appends it again
                    unused = d.unused_data
                    tail = d.flush()
                    if tail:
                        q.put(tail)
                    if unused:
                        comp = unused
                        d = zlib.decompressobj(47)
                        continue
                    comp = fh.read(size)
                    if not comp:
                        break
                    d = zlib.decompressobj(47)
                    continue
                comp = fh.read(size)
            if not d.eof and fh.tell():
                raise EOFError(f"Compressed file '{f}' ended before the end of stream")
    except BaseException as e:  # noqa: BLE001 - re-raised by the caller
        q.put(e)
    finally:
        q.put(None)


def _count_threaded(f: pathlib.Path) -> Counter[str]:
    """Count with zlib inflate (producer thread) overlapped with counting."""
    # Bounded: the producer waits instead of piling up decompressed data
    q: queue.Queue = queue.Queue(maxsize=4)
    reader = _QueueReader(q)
    th = threading.Thread(target=_inflate_stream, args=(f, q), daemon=True)
    th.start()
    counts = _count_header_lines(io.TextIOWrapper(reader, encoding="latin-1"))
    th.join()
    if reader.error is not None:
        raise reader.error
    return counts


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
    boolean is True only when every member boundary falls on a FASTQ record
    boundary (the decompressed content up to each non-final member's end is
    empty or a whole number of 4-line records); only then can a byte range be split at member boundaries and
    counted line by line. One sequential decompress pass; nothing is retained.
    A misaligned boundary aborts early, so a bad file is rejected in seconds.
    """
    offsets = [0]
    size = f.stat().st_size
    d = zlib.decompressobj(47)
    pending = b""
    comp_read = 0
    last_byte: int | None = None  # last decompressed byte of the current member
    newlines = 0  # decompressed line count so far
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
                newlines += out.count(b"\n")
            leftover = d.unused_data
            if d.eof:
                tail = d.flush()
                if tail:
                    last_byte = tail[-1]
                    newlines += tail.count(b"\n")
                end = comp_read - len(leftover)
                offsets.append(end)
                if end < size and (last_byte not in (None, 10) or newlines % 4):
                    return offsets, False
                d = zlib.decompressobj(47)
                last_byte = None
            pending = leftover
    if offsets[-1] != size:
        raise EOFError(f"Compressed file '{f}' ended before the end of stream")
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
        # The range starts on a record boundary, so every 4th line is a header
        for line in itertools.islice(fh, 0, None, 4):
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


def count_header_indexes(
    f: pathlib.Path, budget: int = 1, cache_dir: pathlib.Path | None = None
) -> Counter[str]:
    """Count the header indexes of a gzipped or plain FASTQ file.

    The first header is checked, as the counting pass trusts every header to
    end with ``:<index>``. ``budget`` and ``cache_dir`` select the parallel
    strategies of ``_count_for_file`` (gzipped files only).

    Raises:
        ValueError: If the first record is not a FASTQ header with an index
            (``@<name> <read>:<filter>:<control>:<index>``).
    """
    with _open_text(f) as fh:
        first = fh.readline().rstrip()
    if not first:
        return Counter()
    name, _, comment = first.partition(" ")
    if not name.startswith("@") or ":" not in comment:
        raise ValueError(f"Invalid read name format in {f.name}: {first}")
    if not _is_gzip(f):
        return _count_serial(f)
    return _count_for_file(f, budget, cache_dir)
