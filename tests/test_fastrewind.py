"""Tests for the fastrewind module."""

import argparse

import pytest

from biomate.fastrewind.fastrewind import (
    clean_directory,
    parse_sequence_mask,
    validate_args,
)


# ---------------------------------------------------------------------------
# parse_sequence_mask
# ---------------------------------------------------------------------------


class TestParseSequenceMask:
    def test_r1_only(self):
        """A single-segment Y mask sets R1 and leaves all other fields at 0."""
        result = parse_sequence_mask("Y151")
        assert result["R1"] == 151
        assert result["R2"] == 0

    def test_r1_and_r2(self):
        """A four-segment mask sets both R1 and R2 correctly."""
        # A 2-segment mask "Y151;Y151" maps both to index-1 (R1) due to i_substring//2+1.
        # R2 requires 4 segments; use a full 4-part mask to verify R2.
        result = parse_sequence_mask("Y151;I0;I0;Y151")
        assert result["R1"] == 151
        assert result["R2"] == 151
        assert result["I1"] == 0
        assert result["I2"] == 0

    def test_dual_index(self):
        """A four-segment mask with two I segments sets I1, I2, R1 and R2."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert result["R1"] == 151
        assert result["I1"] == 10
        assert result["I2"] == 10
        assert result["R2"] == 151

    def test_single_index_with_trailing_trim(self):
        """An I8N2 segment sets I1=8 and the trailing-trim I1A=2."""
        result = parse_sequence_mask("Y151;I8N2;Y151")
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_prefix_and_suffix_around_index(self):
        """N before I sets I1B; N after I sets I1A."""
        result = parse_sequence_mask("Y151;N1I8N2;Y151")
        assert result["I1B"] == 1
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_umi_first_position(self):
        """A UMI prefix (U8) before a read segment is parsed into U1."""
        result = parse_sequence_mask("U8Y143;I10;I10;Y151")
        assert result["U1"] == 8
        assert result["R1"] == 143

    def test_too_many_sections_raises(self):
        """A mask with more than 4 semicolon-separated sections raises ValueError."""
        with pytest.raises(ValueError, match="OverrideCycles"):
            parse_sequence_mask("Y151;I10;I10;Y151;Y50")

    def test_all_zero_defaults(self):
        """Unused mask fields default to 0 in the returned dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert result["U1"] == 0
        assert result["U2"] == 0
        assert result["R1B"] == 0
        assert result["R1A"] == 0
        assert result["I1B"] == 0
        assert result["I2B"] == 0

    def test_return_type_is_dict(self):
        """parse_sequence_mask returns a dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        assert isinstance(result, dict)


# ---------------------------------------------------------------------------
# clean_directory
# ---------------------------------------------------------------------------


class TestCleanDirectory:
    def test_removes_a_file(self, tmp_path):
        """A plain file is deleted."""
        f = tmp_path / "test.txt"
        f.write_text("data")
        clean_directory(f)
        assert not f.exists()

    def test_removes_empty_directory(self, tmp_path):
        """An empty directory is deleted."""
        d = tmp_path / "emptydir"
        d.mkdir()
        clean_directory(d)
        assert not d.exists()

    def test_removes_directory_with_files(self, tmp_path):
        """A directory containing files is deleted recursively."""
        d = tmp_path / "mydir"
        d.mkdir()
        (d / "a.txt").write_text("hello")
        (d / "b.txt").write_text("world")
        clean_directory(d)
        assert not d.exists()

    def test_removes_nested_directories(self, tmp_path):
        """A deeply nested directory tree is deleted completely."""
        root = tmp_path / "root"
        root.mkdir()
        nested = root / "a" / "b" / "c"
        nested.mkdir(parents=True)
        (nested / "file.bin").write_bytes(b"\x00\x01")
        clean_directory(root)
        assert not root.exists()


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _make_args(self, tmp_path, force=False):
        input_path = tmp_path / "flowcell"
        input_path.mkdir()
        # Add a dummy FASTQ file so the glob check passes
        (input_path / "sample.fastq.gz").write_bytes(b"")
        # Add a sample sheet
        (input_path / "SampleSheet.csv").write_text("[BCLConvert_Data]\n")
        return argparse.Namespace(
            input_path=input_path,
            output_path=tmp_path / "output",
            sample_sheet=None,
            threads=0,
            force=force,
        )

    def test_valid_args_returned(self, tmp_path):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._make_args(tmp_path)
        result = validate_args(args)
        assert result is args

    def test_nonexistent_input_path_raises(self, tmp_path):
        """A non-existent input path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.input_path = tmp_path / "missing"
        with pytest.raises(argparse.ArgumentTypeError, match="[Vv]alid directory"):
            validate_args(args)

    def test_input_path_is_file_raises(self, tmp_path):
        """A plain file passed as input path raises ArgumentTypeError."""
        f = tmp_path / "not_a_dir.txt"
        f.write_text("")
        args = self._make_args(tmp_path)
        args.input_path = f
        with pytest.raises(argparse.ArgumentTypeError, match="[Vv]alid directory"):
            validate_args(args)

    def test_no_fastq_files_raises(self, tmp_path):
        """A directory with no FASTQ files raises ArgumentTypeError."""
        empty_dir = tmp_path / "empty"
        empty_dir.mkdir()
        (empty_dir / "SampleSheet.csv").write_text("")
        args = self._make_args(tmp_path)
        args.input_path = empty_dir
        with pytest.raises(argparse.ArgumentTypeError, match="[Nn]o FASTQ"):
            validate_args(args)

    def test_missing_sample_sheet_raises(self, tmp_path):
        """Absence of both SampleSheet.csv and --sample-sheet raises ArgumentTypeError."""
        input_path = tmp_path / "flowcell2"
        input_path.mkdir()
        (input_path / "sample.fastq.gz").write_bytes(b"")
        # No SampleSheet.csv and no --sample-sheet flag
        args = argparse.Namespace(
            input_path=input_path,
            output_path=tmp_path / "output",
            sample_sheet=None,
            threads=0,
            force=False,
        )
        with pytest.raises(argparse.ArgumentTypeError, match="[Ss]ample[Ss]heet"):
            validate_args(args)

    def test_explicit_sample_sheet_nonexistent_raises(self, tmp_path):
        """A non-existent --sample-sheet path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.sample_sheet = tmp_path / "nonexistent.csv"
        with pytest.raises(argparse.ArgumentTypeError, match="[Ss]ample"):
            validate_args(args)
