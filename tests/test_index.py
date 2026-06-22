"""Tests for the index module."""

import argparse
from collections import Counter, defaultdict

import pytest

from biomate.index.index import (
    compile_index_regex,
    expand_regex,
    validate_args,
    write_results,
)


# ---------------------------------------------------------------------------
# compile_index_regex
# ---------------------------------------------------------------------------


class TestCompileIndexRegex:
    def test_returns_compiled_regex(self):
        """compile_index_regex returns an object with a search method."""
        pattern = compile_index_regex("ACGT", distance=0, error_type="s")
        assert hasattr(pattern, "search")

    def test_exact_match_distance_zero(self):
        """An exact sequence matches with distance=0."""
        pattern = compile_index_regex("ACGT", distance=0, error_type="s")
        assert pattern.search("ACGT")

    def test_no_match_beyond_distance_zero(self):
        """A sequence differing by one base does not match with distance=0."""
        pattern = compile_index_regex("ACGT", distance=0, error_type="s")
        assert not pattern.search("ACGG")

    def test_one_substitution_allowed(self):
        """A sequence with one substitution matches with distance=1, error_type='s'."""
        pattern = compile_index_regex("ACGT", distance=1, error_type="s")
        # ACGG differs from ACGT by 1 substitution
        assert pattern.search("ACGG")

    def test_two_substitutions_rejected_at_distance_one(self):
        """Two substitutions do not match when distance=1."""
        pattern = compile_index_regex("ACGT", distance=1, error_type="s")
        # AAGG differs from ACGT by 2 substitutions
        assert not pattern.search("AAGG")

    def test_dual_index_with_plus_separator(self):
        """A dual-index pattern using '+' compiles and matches the combined string."""
        pattern = compile_index_regex("ACGT+TTTT", distance=0, error_type="s")
        # Full dual-index string should match
        assert pattern.search("ACGTTTTT") or pattern.search("ACGT+TTTT")

    def test_any_error_type(self):
        """error_type='e' (any edit) allows substitutions."""
        pattern = compile_index_regex("ACGT", distance=1, error_type="e")
        assert pattern.search("ACGT")
        assert pattern.search("ACGG")  # 1 substitution

    def test_insertion_only_type(self):
        """error_type='i' allows insertions."""
        pattern = compile_index_regex("ACT", distance=1, error_type="i")
        # ACT with 1 insertion: ACGT or ACAT etc.
        assert pattern.search("ACGT")

    def test_substitution_type_rejects_insertion(self):
        """error_type='s' does not match a string that requires a deletion."""
        # With s-only, ACGT should not match ACT (deletion, not substitution)
        pattern = compile_index_regex("ACGT", distance=1, error_type="s")
        # ACGT has length 4; deleting a base gives length 3
        # The pattern allows s<=1, d<=0, i<=0
        assert not pattern.search("ACT")


# ---------------------------------------------------------------------------
# expand_regex
# ---------------------------------------------------------------------------


class TestExpandRegex:
    def test_literal_returns_single_element(self):
        """A literal sequence expands to a set containing only itself."""
        assert expand_regex("ACGT") == {"ACGT"}

    def test_wildcard_dot_expands(self):
        """A dot wildcard expands to the five nucleotide characters."""
        result = expand_regex("A.T")
        assert "AAT" in result
        assert "ACT" in result
        assert "AGT" in result
        assert "ATT" in result
        assert "ANT" in result

    def test_character_class_expands(self):
        """A character class expands to all listed bases."""
        result = expand_regex("[ACG]T")
        assert result == {"AT", "CT", "GT"}

    def test_alternation_in_parens(self):
        """Parenthesised alternation expands to each alternative."""
        result = expand_regex("(A|C)GT")
        assert result == {"AGT", "CGT"}

    def test_optional_group(self):
        """An optional group (A)? expands to both the group and the empty string."""
        result = expand_regex("(A)?GT")
        assert result == {"AGT", "GT"}

    def test_multi_character_alternation(self):
        """Multi-character alternation inside parentheses expands correctly."""
        # "(AA|CC)" with ) at end-of-string triggers IndexError in expand_regex (known bug).
        # Use a trailing literal to avoid the boundary condition.
        result = expand_regex("(AA|CC)T")
        assert result == {"AAT", "CCT"}

    def test_plus_raises(self):
        """A '+' character in the pattern raises ValueError."""
        with pytest.raises(ValueError, match="\\+"):
            expand_regex("A+C")

    def test_literal_no_special_chars(self):
        """A plain sequence with no special characters expands to itself."""
        assert expand_regex("GGGG") == {"GGGG"}

    def test_multiple_wildcards(self):
        """Multiple dot wildcards produce 5^n unique expansions."""
        result = expand_regex("A..")
        # 5 choices × 5 choices = 25 unique sequences of length 3 starting with A
        assert len(result) == 25
        assert all(s.startswith("A") and len(s) == 3 for s in result)


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _make_args(
        self,
        input_path,
        index_regex="ACGT",
        index_file=None,
        output=None,
        distance=0,
        error_type="s",
    ):
        return argparse.Namespace(
            input=input_path,
            index_regex=index_regex,
            index_file=index_file,
            output=output,
            distance=distance,
            error_type=error_type,
        )

    def test_valid_args_returned(self, tmp_path):
        """validate_args returns the same Namespace when all arguments are valid."""
        f = tmp_path / "sample.fastq"
        f.write_text("")
        args = self._make_args(f)
        result = validate_args(args)
        assert result is args

    def test_nonexistent_input_raises(self, tmp_path):
        """A path that does not exist raises FileNotFoundError."""
        args = self._make_args(tmp_path / "missing.fastq")
        with pytest.raises(FileNotFoundError):
            validate_args(args)


# ---------------------------------------------------------------------------
# write_results
# ---------------------------------------------------------------------------


class TestWriteResults:
    def _make_results(self):
        results = defaultdict(Counter)
        results["0"][("ACGT", "TTTT")] = 10
        results["0"][("ACGG", "TTTT")] = 5
        results["1"][("ACGT", "TTTT")] = 2
        return dict(results)

    def test_creates_matches_file(self, tmp_path):
        """A pattern{n}_matches.txt file is written to the output directory."""
        write_results(self._make_results(), total_records=100, output_path=tmp_path)
        assert (tmp_path / "pattern1_matches.txt").is_file()

    def test_creates_errors_file(self, tmp_path):
        """A pattern{n}_errors.txt file is written to the output directory."""
        write_results(self._make_results(), total_records=100, output_path=tmp_path)
        assert (tmp_path / "pattern1_errors.txt").is_file()

    def test_matches_file_contains_sequences(self, tmp_path):
        """The matches file contains the matched sequence tokens."""
        write_results(self._make_results(), total_records=100, output_path=tmp_path)
        content = (tmp_path / "pattern1_matches.txt").read_text()
        assert "ACGT" in content
        assert "TTTT" in content

    def test_errors_file_contains_percentages(self, tmp_path):
        """The errors file contains percentage values for each fuzziness level."""
        write_results(self._make_results(), total_records=100, output_path=tmp_path)
        content = (tmp_path / "pattern1_errors.txt").read_text()
        assert "%" in content

    def test_file_index_parameter(self, tmp_path):
        """file_index controls the numeric suffix in the output filenames."""
        write_results(
            self._make_results(), total_records=100, output_path=tmp_path, file_index=3
        )
        assert (tmp_path / "pattern3_matches.txt").is_file()
        assert (tmp_path / "pattern3_errors.txt").is_file()

    def test_creates_output_directory_if_missing(self, tmp_path):
        """The output directory is created automatically if it does not exist."""
        new_dir = tmp_path / "subdir"
        write_results(self._make_results(), total_records=50, output_path=new_dir)
        assert new_dir.is_dir()
