"""Differential parity tests: biomate_kernels vs the production code paths.

IndexPattern must be byte-identical to the `regex`-module path on (match?,
span, fuzzy_counts, groups[1:]); scan_n_content must equal the per-record
dnaio loop in `biomate.nspector.scan`. These tests are the definition of
done for the kernels (see docs/todo/index.md and nspector.md); they skip when
the optional biomate_kernels package is not installed.
"""

import gzip
import pathlib
import random

import pytest

biomate_kernels = pytest.importorskip("biomate_kernels")

from biomate.index.index import compile_index_regex  # noqa: E402
from biomate.nspector.scan import scan_file_python  # noqa: E402

M = "ATCGATCGAT"
G8 = "GGGGGGGG"
ACGT = "ACGT"
SEQ_LEN = 150

PATTERNS = [
    (M, 0, "e"),
    (M, 1, "e"),
    (M, 2, "e"),
    (M, 3, "e"),
    (M, 2, "s"),
    (M, 2, "i"),
    (M, 2, "d"),
    (M, 1, "s"),
    (M, 1, "i"),
    (M, 1, "d"),
    (M, 2, "is"),
    (M, 2, "sd"),
    (M, 3, "id"),
    (f"{M}+{G8}", 2, "e"),
    (f"{M}+{G8}", 1, "e"),
    (f"{M}+{G8}", 2, "s"),
]


def _reference(pypat, seq):
    """Reference result from the production regex engine."""
    m = pypat.search(seq)
    if m is None:
        return None
    return (m.span(), tuple(m.fuzzy_counts), m.groups()[1:])


def _kernel(rustpat, seq):
    """Kernel result in the reference shape."""
    got = rustpat.search(seq)
    if got is None:
        return None
    start, end, counts, parts = got
    return ((start, end), counts, tuple(parts))


def _rand_seq(rng, length=SEQ_LEN):
    return "".join(rng.choice(ACGT) for _ in range(length))


def _plant(rng, base, kind, pos):
    """Insert one 'N' edit into the motif region of a random sequence."""
    seq = list(base)
    start = rng.randrange(0, SEQ_LEN - len(M))
    seg = list(M)
    if kind == "s":
        seg[pos] = "N"
    elif kind == "i":
        seg.insert(pos, "N")
    elif kind == "d":
        del seg[pos]
    seq[start : start + len(M)] = seg
    return "".join(seq)


class TestIndexPatternParity:
    """Kernel vs regex on random sequences across all pattern configs."""

    @pytest.mark.parametrize("motif,distance,error_type", PATTERNS)
    def test_random_sequences(self, motif, distance, error_type):
        rng = random.Random(42)
        pypat = compile_index_regex(motif, distance, error_type)
        rustpat = biomate_kernels.IndexPattern(motif, distance, error_type)
        for _ in range(500):
            seq = _rand_seq(rng)
            assert _kernel(rustpat, seq) == _reference(pypat, seq), seq

    @pytest.mark.parametrize(
        "motif,distance,error_type",
        [(M, 2, "e"), (M, 3, "e"), (M, 2, "s"), (M, 2, "i"), (f"{M}+{G8}", 2, "e")],
    )
    @pytest.mark.parametrize("kind", ["s", "i", "d"])
    def test_planted_edits_all_positions(self, motif, distance, error_type, kind):
        """Single planted edit at every position of the motif."""
        rng = random.Random(7)
        base = _rand_seq(rng)
        pypat = compile_index_regex(motif, distance, error_type)
        rustpat = biomate_kernels.IndexPattern(motif, distance, error_type)
        positions = range(len(M) + (1 if kind == "i" else 0))
        for pos in positions:
            seq = _plant(rng, base, kind, pos)
            assert _kernel(rustpat, seq) == _reference(pypat, seq), (kind, pos, seq)

    def test_known_semantic_cases(self):
        """Probe cases that pin down the engine's documented behaviors."""
        pypat = compile_index_regex(M, 2, "e")
        rustpat = biomate_kernels.IndexPattern(M, 2, "e")
        cases = [
            M,  # exact
            "XXX" + M,  # leftmost start beats fewer errors
            "X" + M,  # non-minimal sub+ins split
            M + M,  # repetition tie-break: shorter span
            "N" + M,
            "NN" + M,
            "ATCGAATCGAT",  # planted insertion
            "ATCGACGAT",  # planted deletion
            "",  # never matches empty
            "ACGTACGT",  # no motif
        ]
        for seq in cases:
            assert _kernel(rustpat, seq) == _reference(pypat, seq), seq


class TestScanBatchParity:
    """scan_batch aggregation must equal the per-record regex aggregation."""

    def test_aggregation_and_fallback(self):
        rng = random.Random(99)
        motif, distance, error_type = M, 2, "e"
        pypat = compile_index_regex(motif, distance, error_type)
        rustpat = biomate_kernels.IndexPattern(motif, distance, error_type)

        seqs = [_rand_seq(rng) for _ in range(500)]
        # one non-ASCII sequence is routed to the fallback list (the caller
        # re-runs those through the regex engine)
        seqs[123] = seqs[123][:50] + "ACGTé" + seqs[123][51:]

        aggregated, fallback = rustpat.scan_batch(seqs)
        assert fallback == [123]

        reference = {}
        for i, seq in enumerate(seqs):
            if i in fallback:
                continue
            m = pypat.search(seq)
            if m:
                key = str(sum(m.fuzzy_counts))
                groups = " ".join(m.groups()[1:])
                reference.setdefault(key, {})
                reference[key][groups] = reference[key].get(groups, 0) + 1

        got = {k: dict(v) for k, v in aggregated.items()}
        assert got == reference

    def test_fallback_indices(self):
        rustpat = biomate_kernels.IndexPattern(M, 2, "e")
        pypat = compile_index_regex(M, 2, "e")
        rng = random.Random(1)
        seqs = [_rand_seq(rng) for _ in range(50)]
        seqs[10] = "ACGTé" + "N" * 50  # non-ASCII char
        seqs[49] = "ACGTé" * 30  # non-ASCII char
        aggregated, fallback = rustpat.scan_batch(seqs)
        assert fallback == [10, 49]
        # every other sequence must still be aggregated (match or not)
        expected_matches = sum(
            1
            for i, seq in enumerate(seqs)
            if i not in (10, 49) and pypat.search(seq) is not None
        )
        total = sum(c for groups in aggregated.values() for c in groups.values())
        assert total == expected_matches


# ---------------------------------------------------------------------------
# nspector N-content scan parity
# ---------------------------------------------------------------------------


def _illumina_name(tile: int, x: int, y: int, index: str = "AA+CC") -> str:
    return f"M00150:2:AABCCC:1:{tile:04d}:{x:04d}:{y:04d} 1:N:0:{index}"


def _write_fastq(path: pathlib.Path, records: list[tuple], gz: bool = True):
    """Records are (name, seq) or (name, seq, qual); qualities default to 'I'."""
    data = "".join(
        f"@{r[0]}\n{r[1]}\n+\n{r[2] if len(r) > 2 else 'I' * len(r[1])}\n"
        for r in records
    )
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(data)
    else:
        path.write_text(data)


def _reference_scan(path: pathlib.Path, bin_size: int = 1000, min_qual: int = 20):
    """Reference: the dnaio loop in nspector.scan.scan_file_python()."""
    scan = scan_file_python(path, bin_size=bin_size, min_qual=min_qual)
    return (
        scan.format_error,
        scan.max_len,
        scan.tiles,
        scan.tile_cycle,
        scan.n_per_read,
        scan.spatial,
        scan.n_shape,
    )


def _normalise(res):
    fmt_err, max_len, *counters = res
    return (fmt_err, max_len, *[[tuple(t) for t in c] for c in counters])


def _kernel_scan(
    path: pathlib.Path, threads: int = 1, bin_size: int = 1000, min_qual: int = 20
):
    return _normalise(
        biomate_kernels.scan_n_content([str(path)], threads, bin_size, min_qual)[0]
    )


class TestScanNContentParity:
    """scan_n_content must equal the dnaio reference scan, counter for counter."""

    def test_mixed_records(self, tmp_path):
        records = [
            (_illumina_name(101, 100, 100), "ACGTNACGTNAC"),  # Ns @5,10
            (_illumina_name(102, 101, 101), "ACGTACGT"),  # no Ns
            (_illumina_name(103, 102, 102), "NNNN"),  # all Ns
            (_illumina_name(104, 1103, 2103), "ACGnTN"),  # lowercase n ignored
            (_illumina_name(105, 104, 104), "ACGTACGTNN"),  # Ns @9,10
        ]
        p = tmp_path / "mixed.fastq.gz"
        _write_fastq(p, records)
        assert _kernel_scan(p) == _reference_scan(p)
        assert _kernel_scan(p, bin_size=7) == _reference_scan(p, bin_size=7)

    def test_format_error_stops_file(self, tmp_path):
        records = [
            (_illumina_name(101, 1, 1), "ACGTNACGT"),  # ok (record 0)
            ("only:six:fields:1:2:3", "ACGTN"),  # error (record 1)
            (_illumina_name(102, 2, 2), "ACGTN"),  # never read
        ]
        p = tmp_path / "bad.fastq.gz"
        _write_fastq(p, records)
        got = _kernel_scan(p)
        assert got == _reference_scan(p)
        assert got[0] == 1  # 0-based index of the offending record

    def test_format_error_without_n(self, tmp_path):
        """Every header is parsed (for totals), so a bad name without Ns also stops."""
        records = [
            (_illumina_name(101, 1, 1), "ACGTNACGT"),
            ("short:bad:name", "ACGT"),
            (_illumina_name(102, 2, 2), "ACGTN"),
        ]
        p = tmp_path / "bad_no_n.fastq.gz"
        _write_fastq(p, records)
        got = _kernel_scan(p)
        assert got == _reference_scan(p)
        assert got[0] == 1

    def test_counters(self, tmp_path):
        records = [
            ("M1:1:FC:1:1101:10:10 1:N:0:AA", "ACGT"),
            ("M1:1:FC:1:1101:1500:11 1:N:0:AA", "ACGN"),
            ("M1:1:FC:2:1101:12:12 1:N:0:AA", "ACGTAC"),
            ("M1:1:FC:2:1102:13:2013 1:N:0:AA", "NN"),
        ]
        p = tmp_path / "lanes.fastq.gz"
        _write_fastq(p, records)
        got = _kernel_scan(p)
        assert got == _reference_scan(p)
        fmt_err, max_len, tiles, tile_cycle, n_per_read, spatial, n_shape = got
        assert (fmt_err, max_len) == (None, 6)
        assert tiles == [
            (1, 1101, 2, 1, 1, 0),
            (2, 1101, 1, 0, 0, 0),
            (2, 1102, 1, 1, 2, 0),
        ]
        assert tile_cycle == [
            (1, 1101, 4, 1, 0),
            (2, 1102, 1, 1, 0),
            (2, 1102, 2, 1, 0),
        ]
        assert n_per_read == [(1, 1, 1), (2, 2, 1)]
        assert spatial == [(1, 1101, 1, 0, 1), (2, 1102, 0, 2, 1)]
        assert n_shape == [(1, 1, 0, 0, 1), (2, 2, 0, 1, 1)]

    def test_no_n_records(self, tmp_path):
        p = tmp_path / "clean.fastq.gz"
        _write_fastq(p, [(_illumina_name(1, 2, 3), "ACGT")])
        got = _kernel_scan(p)
        assert got == _reference_scan(p)
        assert got[3:] == ([], [], [], [])

    def test_multi_member_gzip(self, tmp_path):
        """BCL Convert writes concatenated gzip members; all must be read."""
        p = tmp_path / "multi.fastq.gz"
        with open(p, "wb") as fh:
            for i in range(3):
                fh.write(
                    gzip.compress(
                        f"@{_illumina_name(1101 + i, i, i)}\nACGN\n+\nIIII\n".encode()
                    )
                )
        got = _kernel_scan(p)
        assert got == _reference_scan(p)
        assert sum(t[2] for t in got[2]) == 3

    def test_plain_fastq(self, tmp_path):
        records = [
            (_illumina_name(11, 12, 13), "NNACGT"),
            (_illumina_name(14, 15, 16), "ACGT"),
        ]
        p = tmp_path / "plain.fastq"
        _write_fastq(p, records, gz=False)
        assert _kernel_scan(p) == _reference_scan(p)

    def test_coordinate_formats(self, tmp_path):
        records = [
            ("M1:1:FC:1:0101:1000:2000 1:N:0:AA+CC", "ACGTN"),  # 7 fields
            ("M1:1:FC:1:1101:1000:2000 1:N:0:AA+CC", "NN"),  # 5-digit tile
            ("M1:1:FC:1:1101:1000:2000:9 1:N:0:AA+CC", "N"),  # 8 fields
        ]
        p = tmp_path / "coords.fastq.gz"
        _write_fastq(p, records)
        assert _kernel_scan(p) == _reference_scan(p)

    @pytest.mark.parametrize("min_qual", [0, 2, 20, 30, 93])
    def test_low_quality_and_n_shape(self, tmp_path, min_qual):
        """Random sequences and qualities: low-Q and N-run counters match."""
        rng = random.Random(min_qual)
        records = []
        for _ in range(2_000):
            length = rng.choice((5, 12, 25, 151))
            seq = "".join(rng.choice("ACGTNNN") for _ in range(length))
            if rng.random() < 0.05:
                seq = "N" * length
            qual = "".join(chr(33 + rng.randint(0, 41)) for _ in range(length))
            records.append((_illumina_name(rng.choice((1101, 2204)), 1, 1), seq, qual))
        p = tmp_path / "q.fastq.gz"
        _write_fastq(p, records)
        got = _kernel_scan(p, min_qual=min_qual)
        assert got == _reference_scan(p, min_qual=min_qual)
        assert got[6]  # n_shape populated
        if min_qual in (20, 30):
            assert any(row[4] for row in got[3])  # some low-Q cycles

    def test_large_gzip_spans_many_chunks(self, tmp_path):
        """Exercise the decompress-ahead thread across many 1 MiB chunks."""
        rng = random.Random(3)
        records = [
            (
                _illumina_name(
                    rng.choice((1101, 1102)), rng.randint(0, 9999), rng.randint(0, 9999)
                ),
                "".join(rng.choice("ACGTN") for _ in range(100)),
            )
            for _ in range(30_000)
        ]
        p = tmp_path / "large.fastq.gz"
        _write_fastq(p, records)
        assert _kernel_scan(p) == _reference_scan(p)

    def test_multi_file_threads(self, tmp_path):
        paths = []
        for f in range(3):
            p = tmp_path / f"f{f}.fastq.gz"
            _write_fastq(
                p,
                [
                    (_illumina_name(100 + f * 10 + i, i, i * 2), f"ACGTN{i}ACGT")
                    for i in range(20)
                ],
            )
            paths.append(p)
        got = biomate_kernels.scan_n_content([str(p) for p in paths], 3)
        assert [_normalise(g) for g in got] == [_reference_scan(p) for p in paths]

    def test_bad_bin_size_rejected(self, tmp_path):
        p = tmp_path / "x.fastq.gz"
        _write_fastq(p, [(_illumina_name(1, 2, 3), "ACGT")])
        with pytest.raises(ValueError, match="bin_size"):
            biomate_kernels.scan_n_content([str(p)], 1, 0)
