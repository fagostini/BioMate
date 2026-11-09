"""Tests for the nspector module."""

import argparse
import json
import pathlib
import re

import polars
import pytest

from tests.conftest import illumina_header, write_fastq_gz
from biomate.nspector import nspector
from biomate.nspector.nspector import (
    exploded_frame,
    plot_individual_tiles,
    plot_n_per_cycle_and_tile,
    plot_n_per_sequence,
    usable_scans,
    validate_args,
)
from biomate.nspector.report import (
    build_chart,
    select_spatial,
    write_multiqc,
    write_report,
)
from biomate.nspector.scan import (
    FileScan,
    extract_tile_coordinates,
    fastq_stem,
    parse_fastq_name,
    scan_file_python,
    scan_files,
)
from biomate.nspector.tables import build_tables


def _rec(seq: str, tile="1101", x="100", y="200", lane="1"):
    return (illumina_header(tile=tile, x=x, y=y, lane=lane), seq, "I" * len(seq))


@pytest.fixture()
def run_files(tmp_path) -> list[pathlib.Path]:
    """Two samples, two lanes; tile 1102 on lane 1 has extra Ns late in the read."""
    files = []
    for sample in ("SampA", "SampB"):
        for lane in (1, 2):
            records = []
            for i in range(8):
                tile = "1101" if i % 2 == 0 else "1102"
                seq = "ACGTACGT"
                if i == 0:
                    seq = "NCGTACGT"  # N at cycle 1
                if tile == "1102" and lane == 1:
                    seq = "ACGTACNN"  # Ns at cycles 7 and 8
                records.append(
                    _rec(
                        seq,
                        tile=tile,
                        x=str(100 + i * 10),
                        y=str(500 + i * 7),
                        lane=str(lane),
                    )
                )
            path = tmp_path / f"{sample}_S1_L00{lane}_R1_001.fastq.gz"
            write_fastq_gz(path, records)
            files.append(path)
    return files


# ---------------------------------------------------------------------------
# scan
# ---------------------------------------------------------------------------


class TestScan:
    """Header parsing, file naming and the Python scan."""

    def test_extract_tile_coordinates_includes_lane(self):
        assert extract_tile_coordinates("M1:1:FC:3:1101:10:20 1:N:0:AA") == (
            3,
            1101,
            10,
            20,
        )

    def test_extract_tile_coordinates_invalid(self):
        assert extract_tile_coordinates("bad:name") is None
        assert extract_tile_coordinates("M1:1:FC:x:1101:10:20") is None

    @pytest.mark.parametrize(
        "name, expected",
        [
            ("Sample1_S3_L002_R2_001.fastq.gz", ("Sample1", "R2")),
            ("My_Sample_S12_L001_I1_001.fq.gz", ("My_Sample", "I1")),
            ("Sample1_S3_R1_001.fastq", ("Sample1", "R1")),
            ("whatever.fastq.gz", ("whatever", "R1")),
        ],
    )
    def test_parse_fastq_name(self, name, expected):
        assert parse_fastq_name(pathlib.Path(name)) == expected

    def test_fastq_stem_strips_double_suffix(self):
        assert fastq_stem(pathlib.Path("a/b_R1_001.fastq.gz")) == "b_R1_001"

    def test_scan_counts_totals_and_one_based_cycles(self, tmp_path):
        path = tmp_path / "S_S1_L001_R1_001.fastq.gz"
        write_fastq_gz(
            path,
            [
                _rec("NACGT", tile="1101"),
                _rec("ACGTN", tile="1101"),
                _rec("ACGT", tile="1101"),
                _rec("ACGTAC", tile="1102"),
            ],
        )
        scan = scan_file_python(path)
        assert scan.format_error is None
        assert scan.cycles == [[1], [5]]
        assert scan.tile_totals == {(1, 1101): 3, (1, 1102): 1}
        assert scan.max_len == 6
        assert scan.total_reads == 4
        assert (scan.sample, scan.read) == ("S", "R1")

    def test_scan_stops_on_bad_header(self, tmp_path):
        path = tmp_path / "bad.fastq.gz"
        write_fastq_gz(path, [_rec("ACGT"), ("bad:name", "ACGT", "IIII")])
        scan = scan_file_python(path)
        assert scan.format_error == 1
        assert usable_scans([scan]) == []

    def test_rust_engine_matches_python(self, run_files):
        pytest.importorskip("biomate_kernels")
        py = scan_files(run_files, "python")
        rs = scan_files(run_files, "rust", threads=2)
        assert py == rs


# ---------------------------------------------------------------------------
# tables
# ---------------------------------------------------------------------------


class TestTables:
    """Tidy summary tables and rates."""

    @pytest.fixture()
    def tables(self, run_files):
        return build_tables(scan_files(run_files), spatial_bins=4)

    def test_sample_lane_rates(self, tables):
        row = tables.sample_lane.filter(
            (polars.col("sample") == "SampA") & (polars.col("lane") == 1)
        ).row(0, named=True)
        # read 0 (tile 1101, N at cycle 1) + four tile-1102 reads
        assert (row["n"], row["total"]) == (5, 8)
        assert row["n_rate"] == pytest.approx(5 / 8)

    def test_tile_summary_rate(self, tables):
        row = tables.tile_summary.filter(
            (polars.col("sample") == "SampA")
            & (polars.col("lane") == 1)
            & (polars.col("tile") == 1101)
        ).row(0, named=True)
        assert (row["n"], row["total"]) == (1, 4)
        assert row["n_rate"] == pytest.approx(0.25)
        assert (row["surface"], row["swath"], row["segment"]) == (1, 1, 1)

    def test_tile_without_n_reads_is_kept(self, tables):
        row = tables.tile_summary.filter(
            (polars.col("sample") == "SampA")
            & (polars.col("lane") == 2)
            & (polars.col("tile") == 1102)
        ).row(0, named=True)
        assert (row["n"], row["total"], row["n_rate"]) == (0, 4, 0.0)

    def test_sample_cycle_is_dense(self, tables):
        # 2 samples × 2 lanes × 1 read × 8 cycles
        assert tables.sample_cycle.height == 2 * 2 * 8
        assert tables.sample_cycle["cycle"].min() == 1
        assert tables.sample_cycle["cycle"].max() == 8

    def test_tile_cycle_sums_samples(self, tables):
        row = tables.tile_cycle.filter(
            (polars.col("lane") == 1)
            & (polars.col("tile") == 1102)
            & (polars.col("cycle") == 8)
        ).row(0, named=True)
        # 4 reads per sample on tile 1102, two samples
        assert (row["n"], row["total"]) == (8, 8)

    def test_n_per_read(self, tables):
        hist = tables.n_per_read.filter(
            (polars.col("sample") == "SampA") & (polars.col("lane") == 1)
        )
        assert dict(zip(hist["n_count"], hist["n_reads"])) == {1: 1, 2: 4}

    def test_spatial_bins_in_range(self, tables):
        assert tables.spatial["x_bin"].max() < 4
        assert tables.spatial["y_bin"].max() < 4
        assert tables.spatial["n"].sum() == tables.tile_summary["n"].sum()

    def test_select_spatial_trims_to_worst_tiles(self, tables):
        trimmed = select_spatial(tables, max_rows=1)
        assert trimmed.height <= 1

    def test_write_parquet(self, tables, tmp_path):
        paths = tables.write_parquet(tmp_path)
        assert {p.name for p in paths} >= {
            "nspector_sample_lane.parquet",
            "nspector_tile_cycle.parquet",
            "nspector_spatial.parquet",
        }
        assert (
            polars.read_parquet(tmp_path / "nspector_tile_summary.parquet").height > 0
        )


# ---------------------------------------------------------------------------
# report & MultiQC
# ---------------------------------------------------------------------------


class TestReport:
    """Interactive HTML report and MultiQC custom content."""

    @pytest.fixture()
    def built(self, run_files):
        scans = scan_files(run_files)
        return scans, build_tables(scans, spatial_bins=4)

    def test_chart_has_linked_params(self, built):
        spec = build_chart(built[1]).to_dict()
        names = {p["name"] for p in spec["params"]}
        assert {"lane_sel", "read_sel", "metric"} <= names
        text = json.dumps(spec)
        for selection in ("pick_sample", "pick_tile", "cycle_brush"):
            assert selection in text

    def test_report_is_self_contained(self, built, tmp_path):
        path = write_report(*built, tmp_path)
        html = path.read_text()
        assert "vegaEmbed" in html
        assert not re.search(r"<script[^>]+src=", html)
        assert "<link" not in html

    def test_multiqc_files(self, built, tmp_path):
        paths = write_multiqc(built[1], tmp_path)
        names = {p.name for p in paths}
        assert names == {
            "nspector_tile_rate_mqc.json",
            "nspector_cycle_rate_mqc.json",
            "nspector_general_stats_mqc.tsv",
        }
        heatmap = json.loads(
            (tmp_path / "multiqc" / "nspector_tile_rate_mqc.json").read_text()
        )
        assert heatmap["plot_type"] == "heatmap"
        assert len(heatmap["data"]) == len(heatmap["ycats"])
        assert all(len(row) == len(heatmap["xcats"]) for row in heatmap["data"])
        line = json.loads(
            (tmp_path / "multiqc" / "nspector_cycle_rate_mqc.json").read_text()
        )
        assert line["plot_type"] == "linegraph"
        assert "SampA L001 R1" in line["data"]
        stats = (tmp_path / "multiqc" / "nspector_general_stats_mqc.tsv").read_text()
        assert "SampA\t" in stats and "# plot_type: 'generalstats'" in stats


# ---------------------------------------------------------------------------
# validate_args
# ---------------------------------------------------------------------------


class TestValidateArgs:
    """Tests for validate_args."""

    def _make_args(self, tmp_path, input_files=None, output=None, spatial_bins=50):
        """Create valid input files and a Namespace pointing at them."""
        if input_files is None:
            f = tmp_path / "sample.fastq.gz"
            write_fastq_gz(f, [(illumina_header(), "ACGT", "####")])
            input_files = [f]
        if output is None:
            output = tmp_path / "output"
        return argparse.Namespace(
            input=input_files, output=output, spatial_bins=spatial_bins
        )

    def test_valid_args_returned(self, tmp_path):
        args = self._make_args(tmp_path)
        assert validate_args(args) is args

    def test_nonexistent_input_raises(self, tmp_path):
        args = self._make_args(tmp_path, input_files=[tmp_path / "missing.fastq.gz"])
        with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
            validate_args(args)

    def test_input_is_directory_raises(self, tmp_path):
        d = tmp_path / "adir"
        d.mkdir()
        args = self._make_args(tmp_path, input_files=[d])
        with pytest.raises(argparse.ArgumentTypeError, match="not a file"):
            validate_args(args)

    def test_output_created_when_missing(self, tmp_path):
        new_out = tmp_path / "brand_new"
        validate_args(self._make_args(tmp_path, output=new_out))
        assert new_out.is_dir()

    def test_output_is_file_raises(self, tmp_path):
        out_file = tmp_path / "output.txt"
        out_file.write_text("data")
        args = self._make_args(tmp_path, output=out_file)
        with pytest.raises(argparse.ArgumentTypeError, match="directory"):
            validate_args(args)

    def test_spatial_bins_must_be_positive(self, tmp_path):
        args = self._make_args(tmp_path, spatial_bins=0)
        with pytest.raises(argparse.ArgumentTypeError, match="spatial-bins"):
            validate_args(args)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------


def _main_args(files, output, **kwargs):
    defaults = dict(
        input=files,
        output=output,
        threads=0,
        engine="python",
        spatial_bins=10,
        png=False,
        no_html=False,
        no_multiqc=False,
        verbose=0,
        quiet=False,
        log_file=None,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class TestMain:
    """End-to-end runs of nspector.main."""

    @pytest.fixture(autouse=True)
    def _no_logging_setup(self, monkeypatch):
        monkeypatch.setattr(nspector, "setup_logging", lambda args: None)

    def test_default_outputs(self, run_files, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        nspector.main(_main_args(run_files, out))
        assert (out / "nspector_report.html").is_file()
        assert (out / "multiqc" / "nspector_tile_rate_mqc.json").is_file()
        assert (out / "nspector_tile_cycle.parquet").is_file()
        assert (out / "SampA_S1_L001_R1_001_n_sequences.parquet").is_file()
        assert not (out / "individual_tiles").exists()
        assert not list(out.glob("*.png"))

    def test_opt_outs(self, run_files, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        nspector.main(_main_args(run_files, out, no_html=True, no_multiqc=True))
        assert not (out / "nspector_report.html").exists()
        assert not (out / "multiqc").exists()

    def test_png_opt_in_uses_stem_prefix(self, run_files, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        nspector.main(
            _main_args(run_files[:2], out, png=True, no_html=True, no_multiqc=True)
        )
        names = {p.name for p in (out / "individual_tiles").iterdir()}
        assert "SampA_S1_L001_R1_001_tile_1101.png" in names
        assert "SampA_S1_L002_R1_001_tile_1101.png" in names
        assert (out / "SampA_S1_L001_R1_001_n_count_distribution.png").is_file()


# ---------------------------------------------------------------------------
# static PNG charts
# ---------------------------------------------------------------------------


class TestPlotFunctions:
    """Tests for the PNG plot functions: file creation and determinism."""

    @pytest.fixture()
    def data(self) -> polars.DataFrame:
        scan = FileScan(
            path=pathlib.Path("x_R1_001.fastq.gz"),
            sample="x",
            read="R1",
            lanes=[1] * 6,
            tiles=[1101, 1101, 1101, 1102, 1102, 1103],
            xs=[1, 1, 2, 5, 5, 9],
            ys=[1, 2, 1, 3, 4, 2],
            cycles=[[1, 2], [3], [5], [2], [4, 6], [1]],
            tile_totals={(1, 1101): 3, (1, 1102): 2, (1, 1103): 1},
            max_len=6,
        )
        return exploded_frame(scan)

    def test_plot_n_per_sequence_creates_file(self, tmp_path, data):
        out = tmp_path / "n_count.png"
        plot_n_per_sequence(data, out)
        assert out.stat().st_size > 0

    def test_plot_n_per_cycle_and_tile_creates_file(self, tmp_path, data):
        out = tmp_path / "n_cycle_tile.png"
        plot_n_per_cycle_and_tile(data, out, max_cycle=6)
        assert out.stat().st_size > 0

    def test_plot_individual_tiles_creates_one_file_per_tile(self, tmp_path, data):
        plot_individual_tiles(data, tmp_path, "x_R1_001")
        assert {p.name for p in (tmp_path / "individual_tiles").iterdir()} == {
            "x_R1_001_tile_1101.png",
            "x_R1_001_tile_1102.png",
            "x_R1_001_tile_1103.png",
        }

    def test_plots_are_deterministic(self, tmp_path, data):
        """Rendering the same data twice yields byte-identical PNGs."""
        runs = [tmp_path / "run_a", tmp_path / "run_b"]
        for run in runs:
            run.mkdir()
            plot_n_per_sequence(data, run / "n_count.png")
            plot_n_per_cycle_and_tile(data, run / "n_cycle_tile.png", max_cycle=6)
            plot_individual_tiles(data, run, "x")

        for name in [
            "n_count.png",
            "n_cycle_tile.png",
            "individual_tiles/x_tile_1101.png",
            "individual_tiles/x_tile_1102.png",
            "individual_tiles/x_tile_1103.png",
        ]:
            assert (runs[0] / name).read_bytes() == (runs[1] / name).read_bytes(), name
