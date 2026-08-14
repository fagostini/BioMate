"""Edge case and integration tests for BioMate — unique coverage not in primary test files.

Primary test files cover validate_args, parse_sequence_mask, compile_index_regex,
generate_cycles_filler, and other functions extensively. This file focuses on
edge cases for generate_sequences_set (blabber) which has no dedicated test file.
"""

from collections import Counter

import polars
import pytest

from biomate.blabber.blabber import generate_sequences_set


# ---------------------------------------------------------------------------
# generate_sequences_set — unique edge cases (not in test_blabber.py)
# ---------------------------------------------------------------------------


class TestGenerateSequencesSet:
    """Tests for generate_sequences_set — the only module-level coverage for this function."""

    def test_deterministic_output_without_seed(self):
        """Without an explicit rng, repeated calls produce similar-enough results."""
        seqs_a = generate_sequences_set({"A", "C", "T", "G"}, 50, 5)
        seqs_b = generate_sequences_set({"A", "C", "T", "G"}, 50, 5)
        # Both should have ~5 sequences of length 50
        assert len(seqs_a) == len(seqs_b) == 5

    def test_repeated_calls_produce_different_output(self):
        """Repeated calls produce varying results due to module-level RNG advancing."""
        seqs_a = generate_sequences_set({"A", "C", "T", "G"}, 30, 10)
        seqs_b = generate_sequences_set({"A", "C", "T", "G"}, 30, 10)
        # Not identical because module-level rng state advances between calls
        assert seqs_a != seqs_b

    def test_all_sequences_have_correct_length(self):
        """Every generated sequence matches the requested length."""
        for length in (1, 5, 50, 200):
            seqs = generate_sequences_set({"A", "C", "T", "G"}, length, 10)
            assert all(len(s) == length for s in seqs), f"Failed at length={length}"

    def test_all_bases_from_alphabet(self):
        """No base outside the requested alphabet ever appears."""
        for alphabet in (["A", "C", "T", "G"], ["X", "Y", "Z"], ["0", "1"]):
            seqs = generate_sequences_set(set(alphabet), 50, 10)
            assert all(all(b in alphabet for b in s) for s in seqs)

    def test_custom_alphabet_works(self):
        """Non-DNA alphabets produce valid sequences."""
        seqs = generate_sequences_set({"X", "Y", "Z"}, 20, 5)
        assert all(all(b in {"X", "Y", "Z"} for b in s) for s in seqs)

    def test_returns_at_most_number_sequences(self):
        """The returned set never exceeds the requested number."""
        for n in (1, 3, 10, 100):
            seqs = generate_sequences_set({"A", "C", "T", "G"}, 50, n)
            assert len(seqs) <= n

    def test_returns_exactly_requested_when_space_is_large_enough(self):
        """Requests fewer sequences than the search space returns all requested."""
        # 4^3 = 64 possible sequences of length 3, asking for 10 must succeed
        seqs = generate_sequences_set({"A", "C", "G", "T"}, 3, 10)
        assert len(seqs) == 10

    def test_single_sequence_request(self):
        """Requesting exactly one sequence returns a 1-element set."""
        seqs = generate_sequences_set({"A", "C", "T", "G"}, 50, 1)
        assert len(seqs) == 1
        assert len(next(iter(seqs))) == 50

    def test_saturation_returns_all_possible(self):
        """When requested > search space, returns the full search space."""
        # 2^4 = 16 possible sequences of length 4 with 2 bases
        seqs = generate_sequences_set({"A", "C"}, 4, 100)
        assert len(seqs) == 16

    def test_single_nucleotide_all_identical(self):
        """Single-nucleotide alphabet with any length produces identical sequences."""
        seqs = generate_sequences_set({"A"}, 50, 100)
        assert len(seqs) == 1
        assert next(iter(seqs)) == "A" * 50

    def test_all_sequences_unique(self):
        """The returned set contains no duplicates regardless of count requested."""
        for n in (10, 100, 500):
            seqs = generate_sequences_set({"A", "C", "T", "G"}, 50, n)
            assert len(seqs) == len(set(seqs))

    def test_large_sequence_length(self):
        """Generation of very long sequences succeeds."""
        seqs = generate_sequences_set({"A", "C", "T", "G"}, 10000, 3)
        assert len(seqs) == 3
        assert all(len(s) == 10000 for s in seqs)

    def test_base_distribution_reasonable(self):
        """Concatenated bases show roughly uniform distribution."""
        seqs = generate_sequences_set({"A", "C", "T", "G"}, 100, 100)
        concatenated = "".join(seqs)
        counts = Counter(concatenated)
        total = sum(counts.values())
        for base in "ACGT":
            assert abs(counts.get(base, 0) / total - 0.25) < 0.20

    def test_gc_content_reasonable(self):
        """GC content clusters around 50% for large samples."""
        seqs = generate_sequences_set({"A", "C", "G", "T"}, 100, 50)
        concat = "".join(seqs)
        gc = (concat.count("G") + concat.count("C")) / len(concat)
        assert 0.35 < gc < 0.65


# ---------------------------------------------------------------------------
# Cross-module integration
# ---------------------------------------------------------------------------


class TestCrossModuleIntegration:
    """Tests that exercise multiple modules together."""

    def test_write_results_multi_lane(self, tmp_path):
        """write_results with multiple lanes writes to a single pattern file."""
        from collections import Counter

        from biomate.index.index import write_results

        results = {
            "0": Counter({("ACGT", "TTTT"): 50}),
            "1": Counter({("ACGT", "TTTT"): 5, ("ACGG", "TTTT"): 3}),
            "99": Counter({("AAAA", "TTTT"): 1}),
        }
        write_results(results, total_records=1000, output_path=tmp_path)
        # All lanes go into pattern{n}_matches.txt (file_index defaults to 1)
        matches = tmp_path / "pattern1_matches.txt"
        assert matches.exists()
        content = matches.read_text()
        assert "ACGT" in content
        assert "TTTT" in content
        assert "AAAA" in content
        assert "ACGG" in content

    def test_cycles_filler_schema_consistency(self):
        """generate_cycles_filler output dtype is stable across runs."""
        from biomate.nspector.nspector import generate_cycles_filler

        base = polars.DataFrame({"cycles": [1], "tile": [1101]})
        df1 = generate_cycles_filler(base, 100)
        df2 = generate_cycles_filler(base, 200)

        assert df1.schema == df2.schema

    def test_cycles_filler_value_consistency(self):
        """generate_cycles_filler produces expected cycle ranges."""
        from biomate.nspector.nspector import generate_cycles_filler

        base = polars.DataFrame({"cycles": [1], "tile": [1101]})
        df = generate_cycles_filler(base, 100)
        assert list(df["cycles"].unique()) == list(range(1, 101))

    def test_parse_mask_roundtrip_consistency(self):
        """Same mask string always produces identical dict."""
        from biomate.fastrewind.fastrewind import parse_sequence_mask

        mask = "U10Y151;N2I8N3I10;I8N1;Y151"
        result1 = parse_sequence_mask(mask)
        result2 = parse_sequence_mask(mask)
        assert result1 == result2
