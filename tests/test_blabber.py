"""Tests for the blabber module."""

import argparse

import pytest

from biomate.blabber.blabber import parse_sequence_mask, validate_args


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
        # I8N2 → I1=8, I1A=2 (the N after setting I1 is a suffix)
        result = parse_sequence_mask("Y151;I8N2;Y151")
        assert result["R1"] == 151
        assert result["I1"] == 8
        assert result["I1A"] == 2
        assert result["R2"] == 151

    def test_prefix_and_suffix_around_index(self):
        """N before I sets I1B; N after I sets I1A."""
        # N1I8N2 → I1B=1, I1=8, I1A=2
        result = parse_sequence_mask("Y151;N1I8N2;Y151")
        assert result["I1B"] == 1
        assert result["I1"] == 8
        assert result["I1A"] == 2

    def test_umi_read1(self):
        """A UMI prefix (U8) before a read segment is parsed into U1."""
        # U8Y143 → U1=8, R1=143
        result = parse_sequence_mask("U8Y143;I10;I10;Y151")
        assert result["U1"] == 8
        assert result["R1"] == 143
        assert result["I1"] == 10
        assert result["R2"] == 151

    def test_too_many_sections_raises(self):
        """A mask with more than 4 semicolon-separated sections raises ValueError."""
        with pytest.raises(ValueError, match="OverrideCycles"):
            parse_sequence_mask("Y151;I10;I10;Y151;Y50")

    def test_with_index1_matching_length(self):
        """When index1 length matches the I1 mask value, index1 replaces the int."""
        result = parse_sequence_mask("Y151;I10;I10;Y151", index1="GAACTGAGCG")
        assert result["I1"] == "GAACTGAGCG"

    def test_with_index1_length_mismatch_raises(self):
        """A mismatch between index1 length and the I1 mask value raises ValueError."""
        with pytest.raises(ValueError, match="Index 1"):
            parse_sequence_mask("Y151;I10;I10;Y151", index1="ACGT")  # len 4 ≠ 10

    def test_with_index2_matching_length(self):
        """When index2 length matches the I2 mask value, index2 replaces the int."""
        result = parse_sequence_mask("Y151;I10;I10;Y151", index2="CGCTCCACGA")
        assert result["I2"] == "CGCTCCACGA"

    def test_with_index2_length_mismatch_raises(self):
        """A mismatch between index2 length and the I2 mask value raises ValueError."""
        with pytest.raises(ValueError, match="Index 2"):
            parse_sequence_mask("Y151;I10;I10;Y151", index2="ACGT")  # len 4 ≠ 10

    def test_all_zero_fields_default(self):
        """Unused mask fields default to 0 in the returned dict."""
        result = parse_sequence_mask("Y151;I10;I10;Y151")
        # Unused fields must default to 0
        assert result["U1"] == 0
        assert result["U2"] == 0
        assert result["R1B"] == 0
        assert result["R1A"] == 0


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _base_args(self, **overrides):
        args = argparse.Namespace(
            alphabet="ACGT",
            seq_length=50,
            seq_number=100,
            seq_mask=None,
            index1=None,
            index2=None,
            format="text",
            sample_sheet=None,
            output=None,
        )
        for k, v in overrides.items():
            setattr(args, k, v)
        return args

    def test_valid_args_returned(self):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._base_args()
        assert validate_args(args) is args

    def test_non_alpha_alphabet_raises(self):
        """A non-alphabetic character in --alphabet raises ArgumentTypeError."""
        args = self._base_args(alphabet="ACG1")
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]etters|[Aa]lpha"):
            validate_args(args)

    def test_seq_length_zero_raises(self):
        """--seq-length 0 raises ArgumentTypeError."""
        args = self._base_args(seq_length=0)
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]ength"):
            validate_args(args)

    def test_seq_length_negative_raises(self):
        """A negative --seq-length raises ArgumentTypeError."""
        args = self._base_args(seq_length=-5)
        with pytest.raises(argparse.ArgumentTypeError, match="[Ll]ength"):
            validate_args(args)

    def test_seq_number_zero_raises(self):
        """--seq-number 0 raises ArgumentTypeError."""
        args = self._base_args(seq_number=0)
        with pytest.raises(argparse.ArgumentTypeError, match="[Nn]umber"):
            validate_args(args)

    def test_invalid_index1_non_alpha_raises(self):
        """A non-alphabetic character in --index1 raises ArgumentTypeError."""
        args = self._base_args(index1="ACGT1")
        with pytest.raises(argparse.ArgumentTypeError, match="Index 1"):
            validate_args(args)

    def test_invalid_index2_non_alpha_raises(self):
        """A non-alphabetic character in --index2 raises ArgumentTypeError."""
        args = self._base_args(index2="123")
        with pytest.raises(argparse.ArgumentTypeError, match="Index 2"):
            validate_args(args)

    def test_invalid_seq_mask_raises(self):
        """An invalid --seq-mask raises ValueError."""
        # 5 sections is invalid
        args = self._base_args(seq_mask="Y151;I10;I10;Y151;Y50")
        with pytest.raises(ValueError, match="OverrideCycles"):
            validate_args(args)

    def test_fastq_format_with_sample_sheet_requires_output(self, tmp_path):
        """FASTQ format with --sample-sheet but no --output raises ArgumentTypeError."""
        ss = tmp_path / "ss.csv"
        ss.write_text("")
        args = self._base_args(format="fastq", sample_sheet=ss, output=None)
        with pytest.raises(argparse.ArgumentTypeError, match="[Oo]utput"):
            validate_args(args)
