"""Tests for the strainer module."""

import argparse
import gzip

import polars
import pytest

from tests.conftest import (
    DUAL_INDEX_SS,
    EMPTY_DATA_SS,
    MISSING_COLUMNS_SS,
    NO_DATA_SECTION_SS,
    illumina_header,
    write_fastq_gz,
)
from biomate.strainer.strainer import (
    extract_indexes_from_sample_sheet,
    extract_indexes_from_undetermined_file,
    search_for_unexpected_indexes,
    validate_args,
)


# ---------------------------------------------------------------------------
# extract_indexes_from_sample_sheet
# ---------------------------------------------------------------------------


class TestExtractIndexesFromSampleSheet:
    def test_dual_index_returns_correct_columns(self, dual_index_sample_sheet):
        """Result DataFrame has the expected column names."""
        df = extract_indexes_from_sample_sheet(dual_index_sample_sheet)
        assert set(df.columns) == {"Lane", "Sample_Project", "index", "index2"}

    def test_dual_index_row_count(self, dual_index_sample_sheet):
        """Dual-index sample sheet with 3 entries yields 3 rows."""
        df = extract_indexes_from_sample_sheet(dual_index_sample_sheet)
        assert len(df) == 3

    def test_dual_index_lane_formatting(self, dual_index_sample_sheet):
        """Lane values are formatted as zero-padded strings, e.g. L001."""
        df = extract_indexes_from_sample_sheet(dual_index_sample_sheet)
        assert set(df["Lane"].to_list()) == {"L001", "L002"}

    def test_dual_index_values(self, dual_index_sample_sheet):
        """Index and project values are correctly extracted from the sheet."""
        df = extract_indexes_from_sample_sheet(dual_index_sample_sheet)
        row = df.filter(polars.col("index") == "GAACTGAGCG").row(0, named=True)
        assert row["index2"] == "CGCTCCACGA"
        assert row["Sample_Project"] == "ProjectA"

    def test_single_index_index2_is_none(self, single_index_sample_sheet):
        """Rows from a single-index sheet have index2 as null."""
        df = extract_indexes_from_sample_sheet(single_index_sample_sheet)
        assert df["index2"].is_null().all()

    def test_no_data_section_raises(self, tmp_path):
        """A sheet without a [Data] or [BCLConvert_Data] section raises ValueError."""
        p = tmp_path / "ss.csv"
        p.write_text(NO_DATA_SECTION_SS)
        with pytest.raises(ValueError, match="\\[Data\\]|\\[BCLConvert_Data\\]"):
            extract_indexes_from_sample_sheet(p)

    def test_missing_required_columns_raises(self, tmp_path):
        """A sheet whose header is missing a required column raises ValueError."""
        p = tmp_path / "ss.csv"
        p.write_text(MISSING_COLUMNS_SS)
        with pytest.raises(ValueError, match="Lane"):
            extract_indexes_from_sample_sheet(p)

    def test_empty_data_section_raises(self, tmp_path):
        """A sheet with a header but no data rows raises ValueError."""
        p = tmp_path / "ss.csv"
        p.write_text(EMPTY_DATA_SS)
        with pytest.raises(ValueError, match="No data"):
            extract_indexes_from_sample_sheet(p)

    def test_bclconvert_data_section_detected(self, dual_index_sample_sheet):
        """[BCLConvert_Data] header is accepted."""
        df = extract_indexes_from_sample_sheet(dual_index_sample_sheet)
        assert len(df) > 0

    def test_data_section_v1_also_detected(self, tmp_path):
        """Legacy [Data] section header is parsed in addition to [BCLConvert_Data]."""
        content = DUAL_INDEX_SS.replace("[BCLConvert_Data]", "[Data]")
        p = tmp_path / "ss.csv"
        p.write_text(content)
        df = extract_indexes_from_sample_sheet(p)
        assert len(df) == 3

    def test_lane_single_digit_padded_correctly(self, tmp_path):
        """Lane '1' must become 'L001', lane '12' must become 'L012'."""
        content = (
            "[BCLConvert_Data]\n"
            "Lane,Sample_ID,Sample_Name,index,Sample_Project,OverrideCycles\n"
            "1,S1,S1,ACGT,Proj,Y151;I4;Y151\n"
            "12,S2,S2,TTTT,Proj,Y151;I4;Y151\n"
        )
        p = tmp_path / "ss.csv"
        p.write_text(content)
        df = extract_indexes_from_sample_sheet(p)
        assert "L001" in df["Lane"].to_list()
        assert "L012" in df["Lane"].to_list()


# ---------------------------------------------------------------------------
# extract_indexes_from_undetermined_file
# ---------------------------------------------------------------------------


class TestExtractIndexesFromUndetermined:
    def test_returns_correct_structure(self, undetermined_fastq_l001):
        """Each returned dict has Lane, Index and Count keys."""
        result = extract_indexes_from_undetermined_file(undetermined_fastq_l001)
        assert isinstance(result, list)
        assert all(isinstance(r, dict) for r in result)
        assert all({"Lane", "Index", "Count"} <= r.keys() for r in result)

    def test_lane_extracted_from_filename(self, undetermined_fastq_l001):
        """Lane is extracted from the filename's third underscore-separated token."""
        result = extract_indexes_from_undetermined_file(undetermined_fastq_l001)
        assert all(r["Lane"] == "L001" for r in result)

    def test_counts_are_correct(self, undetermined_fastq_l001):
        """Index occurrence counts reflect the actual frequency in the FASTQ file."""
        result = extract_indexes_from_undetermined_file(undetermined_fastq_l001)
        index_map = {r["Index"]: r["Count"] for r in result}
        # Two reads have the same index → count 2
        assert index_map.get("ACGTACGT+TGCATGCA") == 2

    def test_gggg_indexes_filtered(self, tmp_path):
        """Reads whose index contains 'GGGG' must be excluded."""
        path = tmp_path / "Undetermined_S0_L002_R1_001.fastq.gz"
        write_fastq_gz(
            path,
            [
                (illumina_header(lane="2", index="GGGGACGT+TTTTTTTT"), "ACGT", "####"),
                (illumina_header(lane="2", index="ACGT+TTTT"), "ACGT", "####"),
            ],
        )
        result = extract_indexes_from_undetermined_file(path)
        indexes = [r["Index"] for r in result]
        assert "GGGGACGT+TTTTTTTT" not in indexes
        assert "ACGT+TTTT" in indexes

    def test_invalid_lane_format_raises(self, tmp_path):
        """File whose 3rd underscore-token is not L### must raise ValueError."""
        path = tmp_path / "Undetermined_S0_BADLANE_R1_001.fastq.gz"
        write_fastq_gz(
            path,
            [
                (illumina_header(lane="1", index="ACGT"), "ACGT", "####"),
            ],
        )
        with pytest.raises(ValueError, match="lane"):
            extract_indexes_from_undetermined_file(path)

    def test_malformed_header_raises(self, tmp_path):
        """A header with no space (missing comment section) must raise ValueError."""
        path = tmp_path / "Undetermined_S0_L003_R1_001.fastq.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("@VH00001:1:AABCCC:1:1101:10000:20000\nACGT\n+\n####\n")
        with pytest.raises(ValueError, match="index"):
            extract_indexes_from_undetermined_file(path)

    def test_at_most_1000_entries_returned(self, tmp_path):
        """Result must be capped at 1000 most-common indexes."""
        path = tmp_path / "Undetermined_S0_L004_R1_001.fastq.gz"
        records = [
            (illumina_header(lane="4", index=f"ACGT{i:04d}+TTTT"), "ACGT", "####")
            for i in range(1500)
        ]
        write_fastq_gz(path, records)
        result = extract_indexes_from_undetermined_file(path)
        assert len(result) <= 1000


# ---------------------------------------------------------------------------
# search_for_unexpected_indexes
# ---------------------------------------------------------------------------


def _make_ss_df(rows: list[dict]) -> polars.DataFrame:
    """Build a sample-sheet DataFrame in the format produced by extract_indexes_from_sample_sheet."""
    return polars.from_dicts(rows)


def _make_und_df(rows: list[dict]) -> polars.DataFrame:
    return polars.from_dicts(rows)


class TestSearchForUnexpectedIndexes:
    def test_finds_unexpected_match(self):
        """An index from lane A found in undetermined reads of lane B is returned."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "ProjA",
                    "index": "ACGTACGT",
                    "index2": "TTTTTTTT",
                },
            ]
        )
        und_df = _make_und_df(
            [
                {"Lane": "L002", "Index": "ACGTACGT+TTTTTTTT", "Count": 50},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        assert len(result) == 1
        assert result["Lane"][0] == "L002"
        assert result["Lane_Project"][0] == "L001"

    def test_no_match_returns_empty_typed_df(self):
        """No matches yield an empty typed DataFrame with the correct schema."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "ProjA",
                    "index": "AAAAAAAA",
                    "index2": "TTTTTTTT",
                },
            ]
        )
        und_df = _make_und_df(
            [
                {"Lane": "L002", "Index": "CCCCCCCC+GGGGGGGG", "Count": 10},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        assert result.is_empty()
        assert "Lane" in result.columns
        assert result.schema["Count"] == polars.Int64

    def test_same_lane_not_flagged(self):
        """Undetermined reads from the same lane as the sample must not be flagged."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "ProjA",
                    "index": "ACGTACGT",
                    "index2": "TTTTTTTT",
                },
            ]
        )
        und_df = _make_und_df(
            [
                {"Lane": "L001", "Index": "ACGTACGT+TTTTTTTT", "Count": 100},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        assert result.is_empty()

    def test_one_mismatch_in_index1_still_flagged(self):
        """A single substitution in index1 must still be matched (s<=1)."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "ProjA",
                    "index": "ACGTACGT",
                    "index2": "TTTTTTTT",
                },
            ]
        )
        # ACGAACGT differs from ACGTACGT by 1 substitution
        und_df = _make_und_df(
            [
                {"Lane": "L002", "Index": "ACGAACGT+TTTTTTTT", "Count": 5},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        assert len(result) == 1

    def test_single_index_sample_no_index2(self):
        """Samples with index2=None must still match on index1 alone."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "ProjA",
                    "index": "ACGTACGT",
                    "index2": None,
                },
            ]
        )
        und_df = _make_und_df(
            [
                {"Lane": "L002", "Index": "ACGTACGT", "Count": 20},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        assert len(result) == 1

    def test_empty_schema_columns(self):
        """Empty result must have the same column names as a non-empty one."""
        ss_df = _make_ss_df(
            [
                {
                    "Lane": "L001",
                    "Sample_Project": "P",
                    "index": "AAAA",
                    "index2": "TTTT",
                },
            ]
        )
        und_df = _make_und_df(
            [
                {"Lane": "L002", "Index": "CCCC+GGGG", "Count": 1},
            ]
        )
        result = search_for_unexpected_indexes(ss_df, und_df)
        expected_cols = {
            "Lane",
            "Index",
            "Count",
            "Sample_Project",
            "Lane_Project",
            "index1",
            "index2",
        }
        assert set(result.columns) == expected_cols


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _make_args(self, tmp_path, threads=1, sample_sheet=None):
        fastq = tmp_path / "Undetermined_S0_L001_R1_001.fastq.gz"
        write_fastq_gz(fastq, [(illumina_header(), "ACGT", "####")])
        (tmp_path / "input").mkdir()
        fastq.rename(tmp_path / "input" / "Undetermined_S0_L001_R1_001.fastq.gz")
        (tmp_path / "input" / "SampleSheet.csv").write_text(DUAL_INDEX_SS)
        (tmp_path / "output").mkdir()
        return argparse.Namespace(
            input_path=tmp_path / "input",
            output_path=tmp_path / "output",
            sample_sheet=sample_sheet,
            threads=threads,
        )

    def test_valid_args_returned(self, tmp_path):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._make_args(tmp_path)
        result = validate_args(args)
        assert result is args

    def test_invalid_input_path_raises(self, tmp_path):
        """A non-existent input path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.input_path = tmp_path / "nonexistent"
        with pytest.raises(argparse.ArgumentTypeError):
            validate_args(args)

    def test_threads_zero_raises(self, tmp_path):
        """--threads 0 raises ArgumentTypeError."""
        args = self._make_args(tmp_path, threads=0)
        with pytest.raises(argparse.ArgumentTypeError, match="threads"):
            validate_args(args)

    def test_threads_negative_raises(self, tmp_path):
        """A negative --threads value raises ArgumentTypeError."""
        args = self._make_args(tmp_path, threads=-1)
        with pytest.raises(argparse.ArgumentTypeError, match="threads"):
            validate_args(args)

    def test_output_path_created_if_missing(self, tmp_path):
        """Missing output directory is created automatically during validation."""
        args = self._make_args(tmp_path)
        new_out = tmp_path / "new_output"
        args.output_path = new_out
        validate_args(args)
        assert new_out.is_dir()

    def test_explicit_sample_sheet_validated(self, tmp_path):
        """A non-existent explicit sample sheet path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        args.sample_sheet = tmp_path / "nonexistent.csv"
        with pytest.raises(argparse.ArgumentTypeError, match="[Ss]ample"):
            validate_args(args)
