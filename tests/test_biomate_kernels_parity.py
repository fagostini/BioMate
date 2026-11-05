"""Differential parity tests: biomate_kernels.IndexPattern vs the regex engine.

The kernel must be byte-identical to the production `regex`-module path on
(match?, span, fuzzy_counts, groups[1:]). These tests are the definition of
done for the kernel (see docs/todo/index.md); they skip when the optional
biomate_kernels package is not installed.
"""

import random

import pytest

biomate_kernels = pytest.importorskip("biomate_kernels")

from biomate.index.index import compile_index_regex  # noqa: E402

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
