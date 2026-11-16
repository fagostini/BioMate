"""Tests for the index module."""

import argparse
from collections import Counter, defaultdict

import pytest

from biomate.index.index import (
    compile_index_regex,
    group_indexes,
    main,
    validate_args,
    write_results,
)


# ---------------------------------------------------------------------------
# compile_index_regex
# ---------------------------------------------------------------------------


class TestCompileIndexRegex:
    """Tests for compile_index_regex."""

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
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    """Tests for validate_args."""

    def _make_args(
        self,
        input_path,
        index_regex="ACGT",
        index_file=None,
        output=None,
        distance=0,
        error_type="s",
    ):
        """Build an index Namespace with the given overrides."""
        return argparse.Namespace(
            input=[input_path],
            index_regex=index_regex,
            index_file=index_file,
            output=output,
            distance=distance,
            error_type=error_type,
            min_count=1,
            threads=1,
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

    def test_duplicate_input_names_raise(self, tmp_path):
        """Inputs with the same file name would share an output directory."""
        for d in ("a", "b"):
            (tmp_path / d).mkdir()
            (tmp_path / d / "s.fastq").write_text("")
        args = self._make_args(tmp_path / "a" / "s.fastq")
        args.input.append(tmp_path / "b" / "s.fastq")
        with pytest.raises(ValueError):
            validate_args(args)


# ---------------------------------------------------------------------------
# write_results
# ---------------------------------------------------------------------------


class TestWriteResults:
    """Tests for write_results."""

    def _make_results(self):
        """Build a small results dict across two files and two patterns."""
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


# ---------------------------------------------------------------------------
# main — output without --output, engine choice
# ---------------------------------------------------------------------------


def _main_args(input_path, pattern, engine="regex", output=None):
    """Build a Namespace for index main; input_path may be a list."""
    return argparse.Namespace(
        input=input_path if isinstance(input_path, list) else [input_path],
        index_regex=pattern,
        index_file=None,
        output=output,
        distance=1,
        error_type="s",
        engine=engine,
        min_count=1,
        threads=1,
        verbose=False,
        quiet=False,
    )


def _write_fastq(path, sequences, indexes=None):
    """Write sequences as a plain FASTQ file, optionally with header indexes."""
    names = [
        f"r{i}" if indexes is None else f"r{i} 1:N:0:{indexes[i]}"
        for i in range(len(sequences))
    ]
    path.write_text(
        "".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in zip(names, sequences))
    )


def test_main_prints_results_without_output(tmp_path, capsys):
    """Without --output the matches go to stdout."""
    f = tmp_path / "r.fastq"
    _write_fastq(f, ["AAACGTACGTAA"])
    main(_main_args(f, "ACGTACGT"))
    assert "0 ACGTACGT 1 100.00%" in capsys.readouterr().out


def test_rust_engine_is_skipped_for_regex_syntax(tmp_path, capsys):
    """A pattern with regex syntax gives the regex-engine result on both engines."""
    pytest.importorskip("biomate_kernels")
    f = tmp_path / "r.fastq"
    _write_fastq(f, ["AAACGTACGTAA"])
    main(_main_args(f, "ACG[TA]ACGT", engine="regex"))
    expected = capsys.readouterr().out
    main(_main_args(f, "ACG[TA]ACGT", engine="rust"))
    assert expected and capsys.readouterr().out == expected


# ---------------------------------------------------------------------------
# index discovery and multiple inputs
# ---------------------------------------------------------------------------


class TestGroupIndexes:
    """Tests for group_indexes."""

    def test_one_mismatch_joins_most_frequent_seed(self):
        """A less frequent index within the distance joins the seed."""
        groups = group_indexes(Counter({"ACGT+TTTT": 10, "ACGA+TTTT": 2}), 1)
        assert groups == {"ACGT+TTTT": [(0, "ACGT+TTTT", 10), (1, "ACGA+TTTT", 2)]}

    def test_beyond_distance_or_other_length_seeds_new_group(self):
        """Indexes too far away or of another length start their own group."""
        groups = group_indexes(Counter({"ACGT": 10, "AAAA": 5, "ACGTA": 3}), 1)
        assert list(groups) == ["ACGT", "AAAA", "ACGTA"]

    def test_distance_zero_is_exact_count(self):
        """With distance 0 every distinct index is its own group."""
        groups = group_indexes(Counter({"ACGT": 3, "ACGA": 1}), 0)
        assert groups == {"ACGT": [(0, "ACGT", 3)], "ACGA": [(0, "ACGA", 1)]}


def test_group_indexes_matches_brute_force():
    """Many seeds switch to the variant lookup; results match a plain scan."""
    import random

    rng = random.Random(0)
    counts = Counter(
        {
            "".join(rng.choice("ACG") for _ in range(5)): rng.randint(1, 50)
            for _ in range(150)
        }
    )
    expected: dict = {}
    for index, count in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        for seed, members in expected.items():
            dist = sum(a != b for a, b in zip(seed, index))
            if dist <= 1:
                members.append((dist, index, count))
                break
        else:
            expected[index] = [(0, index, count)]
    assert len(expected) > 20  # more seeds than variants of a 5-mer
    assert group_indexes(counts, 1) == expected


def test_main_discovers_header_indexes(tmp_path, capsys):
    """Without a pattern the header indexes are counted and grouped."""
    f = tmp_path / "r.fastq"
    idx = ["ACGT+TTTT"] * 3 + ["ACGA+TTTT", "GGCC+AATT"]
    _write_fastq(f, ["AAAA"] * 5, idx)
    main(_main_args(f, None))
    assert capsys.readouterr().out.splitlines() == [
        "ACGT+TTTT 0 ACGT+TTTT 3 60.00%",
        "ACGT+TTTT 1 ACGA+TTTT 1 20.00%",
        "GGCC+AATT 0 GGCC+AATT 1 20.00%",
    ]


def test_main_keeps_results_per_input_file(tmp_path, capsys):
    """Several inputs give one block / output subdirectory per file."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    _write_fastq(a, ["AAAA"], ["ACGT"])
    _write_fastq(b, ["AAAA"], ["GGCC"])
    main(_main_args([a, b], None))
    assert capsys.readouterr().out.splitlines() == [
        f"# {a}",
        "ACGT 0 ACGT 1 100.00%",
        f"# {b}",
        "GGCC 0 GGCC 1 100.00%",
    ]

    out = tmp_path / "out"
    main(_main_args([a, b], None, output=out))
    assert (out / "a" / "indexes_matches.txt").read_text() == "ACGT 0 ACGT 1 100.00%\n"
    assert (out / "b" / "indexes_groups.txt").read_text() == "GGCC 1 100.00%\n"


def test_main_min_count_drops_rare_indexes(tmp_path, capsys):
    """Indexes below --min-count are left out and reported as one total."""
    f = tmp_path / "r.fastq"
    idx = ["ACGT+TTTT"] * 3 + ["ACGA+TTTT", "GGCC+AATT"]
    _write_fastq(f, ["AAAA"] * 5, idx)
    args = _main_args(f, None)
    args.min_count = 2
    main(args)
    assert capsys.readouterr().out.splitlines() == [
        "ACGT+TTTT 0 ACGT+TTTT 3 60.00%",
        "# 2 reads in 2 indexes below --min-count 2 (40.00%)",
    ]


def test_main_threads_give_same_files(tmp_path):
    """Parallel inputs write the same per-file results as a serial run."""
    a, b = tmp_path / "a.fastq", tmp_path / "b.fastq"
    _write_fastq(a, ["AAAA"] * 2, ["ACGT", "ACGA"])
    _write_fastq(b, ["AAAA"], ["GGCC"])
    for threads in (1, 2):
        args = _main_args([a, b], None, output=tmp_path / f"t{threads}")
        args.threads = threads
        main(args)
    for name in ("a", "b"):
        for out in ("indexes_matches.txt", "indexes_groups.txt"):
            serial = (tmp_path / "t1" / name / out).read_text()
            assert serial
            assert (tmp_path / "t2" / name / out).read_text() == serial
