"""Tests for the nspector module."""

import argparse

import polars
import pytest

from tests.conftest import illumina_header, write_fastq_gz
from biomate.nspector.nspector import generate_cycles_filler, validate_args


# ---------------------------------------------------------------------------
# generate_cycles_filler
# ---------------------------------------------------------------------------


def _data_df(tiles, x_coords, y_coords, cycles_lists):
    """Build the exploded data DataFrame that generate_cycles_filler expects."""
    return (
        polars.DataFrame(
            {
                "tile": tiles,
                "x": x_coords,
                "y": y_coords,
                "cycles": cycles_lists,
            }
        )
        .explode("cycles")
        .with_columns(
            polars.col("tile").cast(polars.UInt16),
            polars.col("x").cast(polars.UInt32),
            polars.col("y").cast(polars.UInt32),
            polars.col("cycles").cast(polars.UInt16),
        )
    )


class TestGenerateCyclesFiller:
    def test_row_count_is_tiles_times_cycles(self):
        """Filler has exactly (unique tile count) × cycles rows."""
        data = _data_df([1101, 1101, 1102], [0, 0, 0], [0, 0, 0], [[2], [5], [3]])
        result = generate_cycles_filler(data, 5)
        unique_tiles = data["tile"].n_unique()
        assert len(result) == unique_tiles * 5

    def test_n_filler_is_always_zero(self):
        """The n_filler column is 0 for every row."""
        data = _data_df([1101, 1102], [0, 0], [0, 0], [[1], [3]])
        result = generate_cycles_filler(data, 3)
        assert (result["n_filler"] == 0).all()

    def test_cycles_range_is_1_to_max(self):
        """Cycle values span 1 through the supplied max inclusive."""
        data = _data_df([1101], [0], [0], [[4]])
        result = generate_cycles_filler(data, 4)
        assert set(result["cycles"].to_list()) == {1, 2, 3, 4}

    def test_all_unique_tiles_present(self):
        """Every unique tile in the input appears in the filler."""
        data = _data_df(
            [1101, 1101, 1102, 1103], [0] * 4, [0] * 4, [[1], [2], [3], [1]]
        )
        result = generate_cycles_filler(data, 3)
        assert set(result["tile"].to_list()) == {1101, 1102, 1103}

    def test_single_tile_single_cycle(self):
        """A single tile with a single max cycle produces exactly one row."""
        data = _data_df([1101], [0], [0], [[1]])
        result = generate_cycles_filler(data, 1)
        assert len(result) == 1
        assert result["cycles"][0] == 1
        assert result["n_filler"][0] == 0

    def test_columns_present(self):
        """Output columns are exactly tile, cycles and n_filler."""
        data = _data_df([1101], [0], [0], [[2]])
        result = generate_cycles_filler(data, 2)
        assert set(result.columns) == {"tile", "cycles", "n_filler"}


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    def _make_args(self, tmp_path, input_files=None, output=None):
        if input_files is None:
            f = tmp_path / "sample.fastq.gz"
            write_fastq_gz(f, [(illumina_header(), "ACGT", "####")])
            input_files = [f]
        if output is None:
            output = tmp_path / "output"
        return argparse.Namespace(input=input_files, output=output)

    def test_valid_args_returned(self, tmp_path):
        """validate_args returns the same Namespace when all arguments are valid."""
        args = self._make_args(tmp_path)
        result = validate_args(args)
        assert result is args

    def test_nonexistent_input_raises(self, tmp_path):
        """A path that does not exist raises ArgumentTypeError."""
        args = self._make_args(tmp_path, input_files=[tmp_path / "missing.fastq.gz"])
        with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
            validate_args(args)

    def test_input_is_directory_raises(self, tmp_path):
        """Passing a directory as input raises ArgumentTypeError."""
        d = tmp_path / "adir"
        d.mkdir()
        args = self._make_args(tmp_path, input_files=[d])
        with pytest.raises(argparse.ArgumentTypeError, match="not a file"):
            validate_args(args)

    def test_output_created_when_missing(self, tmp_path):
        """Missing output directory is created automatically."""
        new_out = tmp_path / "brand_new"
        args = self._make_args(tmp_path, output=new_out)
        validate_args(args)
        assert new_out.is_dir()

    def test_existing_output_dir_accepted(self, tmp_path):
        """An existing output directory is accepted without error."""
        out = tmp_path / "existing"
        out.mkdir()
        args = self._make_args(tmp_path, output=out)
        result = validate_args(args)
        assert result is args

    def test_output_is_file_raises(self, tmp_path):
        """A plain file used as the output path raises ArgumentTypeError."""
        out_file = tmp_path / "output.txt"
        out_file.write_text("data")
        args = self._make_args(tmp_path, output=out_file)
        with pytest.raises(argparse.ArgumentTypeError, match="directory"):
            validate_args(args)

    def test_multiple_inputs_validated(self, tmp_path):
        """All files in a multi-file input list are validated successfully."""
        f1 = tmp_path / "a.fastq.gz"
        f2 = tmp_path / "b.fastq.gz"
        write_fastq_gz(f1, [(illumina_header(), "ACGT", "####")])
        write_fastq_gz(f2, [(illumina_header(), "TTTT", "####")])
        args = self._make_args(tmp_path, input_files=[f1, f2])
        result = validate_args(args)
        assert result is args
