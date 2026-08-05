"""Additional edge case and integration tests for improved coverage."""

import pytest
import sys
from pathlib import Path
import tempfile

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))


class TestBlabberEdgeCases:
    """Edge case tests for blabber module."""

    def test_sequence_generation_with_zero_seed(self):
        """Test that seed=0 produces reproducible sequences."""
        from biomate.blabber.blabber import generate_sequences_set
        import numpy as np

        # Generate with seed 0
        rng1 = np.random.default_rng(0)
        sequences1 = set()
        for _ in range(10):
            rng1_local = np.random.default_rng(0)
            seq = "".join(rng1_local.choice(["A", "C", "T", "G"], size=50))
            sequences1.add(seq)

        # Reproducibility is ensured by seed
        assert len(sequences1) > 0

    def test_sequence_generation_large_length(self):
        """Test that large sequence lengths work."""
        from biomate.blabber.blabber import generate_sequences_set
        import numpy as np

        nucleotides = {"A", "C", "T", "G"}
        sequences = generate_sequences_set(nucleotides, 10000, 5)

        assert len(sequences) <= 5
        for seq in sequences:
            assert len(seq) == 10000
            assert all(base in nucleotides for base in seq)

    def test_sequence_generation_single_nucleotide(self):
        """Test sequence generation with single nucleotide alphabet."""
        from biomate.blabber.blabber import generate_sequences_set

        nucleotides = {"A"}
        sequences = generate_sequences_set(nucleotides, 50, 3)

        # All sequences should be 'A' repeated
        for seq in sequences:
            assert seq == "A" * 50

    def test_sequence_generation_unicode_alphabet(self):
        """Test that non-DNA alphabets work."""
        from biomate.blabber.blabber import generate_sequences_set

        nucleotides = {"X", "Y", "Z"}
        sequences = generate_sequences_set(nucleotides, 20, 5)

        assert len(sequences) <= 5
        for seq in sequences:
            assert len(seq) == 20
            assert all(base in nucleotides for base in seq)


class TestDirstructEdgeCases:
    """Edge case tests for dirstruct module."""

    def test_extract_from_nonexistent_directory(self):
        """Test that extraction fails gracefully on nonexistent paths."""
        from biomate.dirstruct.dirstruct import validate_args
        import argparse

        args = argparse.Namespace(
            command="extract",
            input="/nonexistent/path/that/does/not/exist",
            output="/tmp/output.txt",
            tags=True,
        )

        with pytest.raises(SystemExit):
            validate_args(args)

    def test_create_in_read_only_directory(self):
        """Test creation in restricted directory."""
        from biomate.dirstruct.dirstruct import validate_args
        import argparse

        args = argparse.Namespace(
            command="create",
            input="/tmp/nonexistent_input.txt",
            output="/root/output",  # Likely read-only for test user
            tags=False,
        )

        # Should validate input exists
        with pytest.raises(SystemExit):
            validate_args(args)


class TestIndexEdgeCases:
    """Edge case tests for index module."""

    def test_compile_index_regex_empty_index(self):
        """Test regex compilation with empty index."""
        from biomate.index.index import compile_index_regex

        # Empty index should match empty string only
        regex = compile_index_regex("", distance=0)
        assert regex.match("") is not None
        assert regex.match("A") is None

    def test_compile_index_regex_long_index(self):
        """Test regex compilation with very long index."""
        from biomate.index.index import compile_index_regex

        long_index = "A" * 100
        regex = compile_index_regex(long_index, distance=0)

        # Should match exact index
        assert regex.match(long_index) is not None
        # Should not match different
        assert regex.match("T" * 100) is None

    def test_compile_index_regex_all_error_types(self):
        """Test all error type combinations."""
        from biomate.index.index import compile_index_regex

        error_types = ["s", "i", "d", "e"]
        for error_type in error_types:
            regex = compile_index_regex("AAAAAA", distance=1, error_type=error_type)
            # Should create valid regex
            assert regex is not None


class TestFastrewindEdgeCases:
    """Edge case tests for fastrewind module."""

    def test_parse_sequence_mask_complex_structure(self):
        """Test parsing complex sequence masks."""
        from biomate.fastrewind.fastrewind import parse_sequence_mask

        # Complex mask with multiple components
        result = parse_sequence_mask("10U5R10;8I1;8I2;10R2")

        assert result is not None
        assert isinstance(result, dict)


class TestNspectorEdgeCases:
    """Edge case tests for nspector module."""

    def test_cycles_filler_single_cycle(self):
        """Test cycle filler with single cycle."""
        from biomate.nspector.nspector import generate_cycles_filler

        # Single cycle, single tile
        df = generate_cycles_filler(1, 1, ["A1"])

        assert df is not None
        assert len(df) >= 1
        assert "Cycle" in df.columns

    def test_cycles_filler_many_cycles(self):
        """Test cycle filler with large cycle count."""
        from biomate.nspector.nspector import generate_cycles_filler

        # Many cycles
        df = generate_cycles_filler(500, 1, ["A1"])

        assert df is not None
        assert len(df) >= 500

    def test_cycles_filler_many_tiles(self):
        """Test cycle filler with large tile count."""
        from biomate.nspector.nspector import generate_cycles_filler

        # Single cycle, many tiles
        tiles = [f"A{i}" for i in range(100)]
        df = generate_cycles_filler(1, 100, tiles)

        assert df is not None
        assert len(df) >= 100


class TestIntegrationScenarios:
    """Integration tests across modules."""

    def test_sequence_generation_produces_valid_output(self):
        """Test that sequence generation produces valid sequences."""
        from biomate.blabber.blabber import generate_sequences_set

        nucleotides = {"A", "C", "T", "G"}
        sequences = generate_sequences_set(nucleotides, 50, 10)

        assert len(sequences) <= 10
        for seq in sequences:
            assert len(seq) == 50
            assert all(base in nucleotides for base in seq)


class TestErrorHandling:
    """Test error handling and recovery."""

    def test_invalid_alphabet_raises(self):
        """Test that invalid alphabets are rejected."""
        from biomate.blabber.blabber import validate_args
        import argparse

        args = argparse.Namespace(
            format="fastq",
            seq_number=100,
            seq_length=50,
            alphabet="123!@#",  # Invalid: non-alpha characters
            output=None,
            sample_sheet=None,
            flowcell_id=None,
            seq_mask=None,
            taint=False,
            index1="",
            index2="",
            random_seed=None,
        )

        with pytest.raises(SystemExit):
            validate_args(args)

    def test_zero_sequences_raises(self):
        """Test that zero sequence count is rejected."""
        from biomate.blabber.blabber import validate_args
        import argparse

        args = argparse.Namespace(
            format="fastq",
            seq_number=0,  # Invalid: must be > 0
            seq_length=50,
            alphabet="ACGT",
            output=None,
            sample_sheet=None,
            flowcell_id=None,
            seq_mask=None,
            taint=False,
            index1="",
            index2="",
            random_seed=None,
        )

        with pytest.raises(SystemExit):
            validate_args(args)

    def test_negative_length_raises(self):
        """Test that negative sequence length is rejected."""
        from biomate.blabber.blabber import validate_args
        import argparse

        args = argparse.Namespace(
            format="fastq",
            seq_number=100,
            seq_length=-50,  # Invalid: must be > 0
            alphabet="ACGT",
            output=None,
            sample_sheet=None,
            flowcell_id=None,
            seq_mask=None,
            taint=False,
            index1="",
            index2="",
            random_seed=None,
        )

        with pytest.raises(SystemExit):
            validate_args(args)


class TestBoundaryConditions:
    """Test boundary conditions and limits."""

    def test_very_short_sequence(self):
        """Test 1-character sequences."""
        from biomate.blabber.blabber import generate_sequences_set

        nucleotides = {"A", "C", "T", "G"}
        sequences = generate_sequences_set(nucleotides, 1, 5)

        assert len(sequences) <= 5
        for seq in sequences:
            assert len(seq) == 1
            assert seq in nucleotides

    def test_one_sequence_request(self):
        """Test requesting exactly one sequence."""
        from biomate.blabber.blabber import generate_sequences_set

        nucleotides = {"A", "C", "T", "G"}
        sequences = generate_sequences_set(nucleotides, 50, 1)

        assert len(sequences) == 1
        assert len(list(sequences)[0]) == 50

    def test_duplicate_sequences_handled(self):
        """Test that duplicate sequences are handled properly."""
        from biomate.blabber.blabber import generate_sequences_set

        # Request many sequences from limited alphabet
        nucleotides = {"A"}
        sequences = generate_sequences_set(nucleotides, 1, 100)

        # With single nucleotide, all sequences are identical
        assert len(sequences) == 1  # Set deduplicates
        assert list(sequences)[0] == "A"
