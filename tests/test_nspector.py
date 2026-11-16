"""Tests for the nspector module."""

import argparse
import json
import pathlib
import random
import re
import sys

import polars
import pytest
import vl_convert

from tests.conftest import illumina_header, write_fastq_gz
from biomate.nspector import nspector
from biomate.nspector import scan as nspector_scan
from biomate.nspector.nspector import (
    plot_individual_tiles,
    plot_n_per_cycle_and_tile,
    plot_n_per_sequence,
    usable_scans,
    validate_args,
)
from biomate.nspector.report import (
    build_chart,
    project_summary,
    select_spatial,
    window_tile_cycle,
    with_tile_anchors,
    write_multiqc,
    write_report,
)
from biomate.nspector.scan import (
    FileScan,
    LowQCounter,
    discover_fastqs,
    extract_tile_coordinates,
    fastq_stem,
    parse_fastq_name,
    project_of,
    scan_file_python,
    scan_files,
)
from biomate.nspector.flags import build_flags, dense_tile_cycle, sample_status
from biomate.nspector.runmetrics import (
    INTEROP_SCHEMA,
    join_tile_interop,
    read_tile_metrics,
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
        assert scan.tile_cycle == [(1, 1101, 1, 1, 0), (1, 1101, 5, 1, 0)]
        # last field: bases sequenced (sum of read lengths)
        assert scan.tiles == [(1, 1101, 3, 2, 2, 0, 14), (1, 1102, 1, 0, 0, 0, 6)]
        assert scan.max_len == 6
        assert scan.total_reads == 4
        assert (scan.sample, scan.read) == ("S", "R1")

    def test_scan_low_quality_and_n_shape(self, tmp_path):
        path = tmp_path / "S_S1_L001_R1_001.fastq.gz"
        seq_a = "ACGT" * 5 + "NNN" + "ACGTACGTACGT"  # run of 3 in the middle
        write_fastq_gz(
            path,
            [
                # Q2 ('#') at cycles 1 and 3; Q19 ('4') at cycle 4; N at 2 is not low-Q
                (illumina_header(), "ANGTA", "##I4I"),
                (illumina_header(), "AC", "5I"),  # Q20 is not low-Q
                (illumina_header(), seq_a, "I" * len(seq_a)),
                (illumina_header(), "NNNN", "####"),
            ],
        )
        scan = scan_file_python(path, min_qual=20)
        assert scan.tiles == [(1, 1101, 4, 3, 8, 1, 46)]
        lowq = {c: q for _, _, c, _, q in scan.tile_cycle if q}
        assert lowq == {1: 1, 4: 1}
        assert scan.n_shape == [
            (1, 1, 0, 0, 1),  # ANGTA: first N at cycle 2
            (1, 3, 1, 0, 1),  # middle run of 3
            (1, 4, 0, 1, 1),  # all N
        ]
        assert scan_file_python(path, min_qual=0).tiles[0][5] == 0

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
        return build_tables(scan_files(run_files, bin_size=25), bin_size=25)

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

    def test_spatial_bins_use_coordinate_width(self, tables):
        # fixture x = 100..170, y = 500..549, bin width 25
        assert set(tables.spatial["x_bin"]) <= {4, 5, 6}
        assert set(tables.spatial["y_bin"]) <= {20, 21}
        assert tables.spatial["n_bases"].sum() == tables.tile_summary["n_bases"].sum()
        # default windows of 10 cycles: Ns at cycles 1, 7 and 8 all start at 1
        assert set(tables.spatial["cycle"]) == {1}

    def test_n_base_rate(self, tables):
        row = tables.sample_lane.row(0, named=True)
        assert row["bases"] == 8 * row["total"]
        assert row["n_base_rate"] == row["n_bases"] / row["bases"]

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
        return scans, build_tables(scans, bin_size=25)

    def test_chart_has_linked_params(self, built):
        spec = build_chart(built[1]).to_dict()
        names = {p["name"] for p in spec["params"]}
        assert {"project_sel", "lane_sel", "read_sel", "metric"} <= names
        selects = {p["name"]: p for p in spec["params"] if "bind" in p}
        assert selects["lane_sel"]["value"] == "All"
        assert selects["lane_sel"]["bind"]["options"][0] == "All"
        assert selects["project_sel"]["bind"]["options"][0] == "All"
        assert selects["show_warn"]["value"] is False
        text = json.dumps(spec)
        assert '"type": "symlog"' in text and "No N bases in this selection" in text
        text = json.dumps(spec)
        for selection in ("pick_sample", "pick_tile", "cycle_brush"):
            assert selection in text

    def test_many_samples_use_median_band(self, built, monkeypatch):
        from biomate.nspector import report

        monkeypatch.setattr(report, "MANY_SAMPLES", 1)
        text = json.dumps(build_chart(built[1]).to_dict())
        assert '"op": "median"' in text and "q3" in text

    def test_tile_cycle_windows_when_large(self, built):
        from biomate.nspector.report import window_tile_cycle

        tables = built[1]
        assert tables.tile_cycle.height == 4
        frame, window = window_tile_cycle(tables, max_rows=3)
        assert window > 1 and frame.height == 3
        # 3 lane/tile groups cannot fit in 2 rows: stops at one window per tile
        frame, window = window_tile_cycle(tables, max_rows=2)
        assert window == 8 and frame.height == 3
        same, one = window_tile_cycle(tables)
        assert one == 1 and same.height == tables.tile_cycle.height

    def test_report_above_altair_row_limit(self, tmp_path):
        """More than Altair's default 5000 rows must not raise MaxRowsError."""
        path = tmp_path / "S_S1_L001_R1_001.fastq.gz"
        write_fastq_gz(path, [_rec("NA", x=str(i), y=str(i)) for i in range(6000)])
        scans = scan_files([path], "python", bin_size=1)
        tables = build_tables(scans, bin_size=1)
        assert tables.spatial.height == 6000
        assert write_report(scans, tables, tmp_path).is_file()

    def test_report_uses_canvas(self, built, tmp_path):
        html = write_report(*built, tmp_path).read_text()
        assert '"renderer": "canvas"' in html

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

    def _make_args(
        self, tmp_path, input_files=None, output=None, spatial_bin_size=1000
    ):
        """Create valid input files and a Namespace pointing at them."""
        if input_files is None:
            f = tmp_path / "sample.fastq.gz"
            write_fastq_gz(f, [(illumina_header(), "ACGT", "####")])
            input_files = [f]
        if output is None:
            output = tmp_path / "output"
        return argparse.Namespace(
            input=input_files, output=output, spatial_bin_size=spatial_bin_size
        )

    def test_valid_args_returned(self, tmp_path):
        args = self._make_args(tmp_path)
        assert validate_args(args) is args

    def test_nonexistent_input_raises(self, tmp_path):
        args = self._make_args(tmp_path, input_files=[tmp_path / "missing.fastq.gz"])
        with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
            validate_args(args)

    def test_directory_without_fastq_raises(self, tmp_path):
        d = tmp_path / "adir"
        d.mkdir()
        (d / "notes.txt").write_text("x")
        args = self._make_args(tmp_path, input_files=[d])
        with pytest.raises(argparse.ArgumentTypeError, match="No FASTQ"):
            validate_args(args)

    def test_directory_is_expanded(self, tmp_path):
        d = tmp_path / "run" / "ProjA"
        d.mkdir(parents=True)
        f = d / "S1_S1_L001_R1_001.fastq.gz"
        write_fastq_gz(f, [(illumina_header(), "ACGT", "####")])
        args = validate_args(self._make_args(tmp_path, input_files=[tmp_path / "run"]))
        assert args.input == [f]

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

    def test_rust_engine_without_kernel_raises(self, tmp_path, monkeypatch):
        monkeypatch.setattr(nspector_scan, "biomate_kernels", None)
        fq = tmp_path / "x.fastq.gz"
        fq.touch()
        args = argparse.Namespace(
            input=[fq], output=tmp_path / "out", engine="rust", spatial_bin_size=1000
        )
        with pytest.raises(argparse.ArgumentTypeError, match="biomate_kernels"):
            validate_args(args)

    def test_spatial_bin_size_must_be_positive(self, tmp_path):
        args = self._make_args(tmp_path, spatial_bin_size=0)
        with pytest.raises(argparse.ArgumentTypeError, match="spatial-bin-size"):
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
        spatial_bin_size=25,
        raw=False,
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
        assert not (out / "SampA_S1_L001_R1_001_n_sequences").exists()
        assert not (out / "individual_tiles").exists()
        assert not list(out.glob("*.png"))

    def test_raw_streams_parts(self, run_files, tmp_path):
        out = tmp_path / "out"
        out.mkdir()
        nspector.main(
            _main_args(run_files[:1], out, raw=True, no_html=True, no_multiqc=True)
        )
        raw = polars.read_parquet(
            out / "SampA_S1_L001_R1_001_n_sequences" / "*.parquet"
        )
        # 1 read with an N at cycle 1 plus 4 tile-1102 reads with Ns at 7 and 8
        assert raw.height == 5
        assert sorted(raw["cycles"].to_list()) == [[1], [7, 8], [7, 8], [7, 8], [7, 8]]

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
    def scan(self, tmp_path) -> FileScan:
        path = tmp_path / "x_R1_001.fastq.gz"
        records = [
            _rec("NNACGT", tile="1101", x="1", y="1"),
            _rec("ACNCGT", tile="1101", x="1", y="2"),
            _rec("ACGTNC", tile="1101", x="2", y="1"),
            _rec("ANGTAC", tile="1102", x="5", y="3"),
            _rec("ACGNAN", tile="1102", x="5", y="4"),
            _rec("NCGTAC", tile="1103", x="9", y="2"),
            _rec("ACGTAC", tile="1103", x="9", y="9"),
        ]
        write_fastq_gz(path, records)
        return scan_file_python(path, bin_size=2)

    def test_plot_n_per_sequence_creates_file(self, tmp_path, scan):
        out = tmp_path / "n_count.png"
        plot_n_per_sequence(scan, out)
        assert out.stat().st_size > 0

    def test_plot_n_per_cycle_and_tile_creates_file(self, tmp_path, scan):
        out = tmp_path / "n_cycle_tile.png"
        plot_n_per_cycle_and_tile(scan, out)
        assert out.stat().st_size > 0

    def test_plot_individual_tiles_creates_one_file_per_tile(self, tmp_path, scan):
        plot_individual_tiles(scan, tmp_path, 2)
        assert {p.name for p in (tmp_path / "individual_tiles").iterdir()} == {
            "x_R1_001_tile_1101.png",
            "x_R1_001_tile_1102.png",
            "x_R1_001_tile_1103.png",
        }

    def test_plots_are_deterministic(self, tmp_path, scan):
        """Rendering the same data twice yields byte-identical PNGs."""
        runs = [tmp_path / "run_a", tmp_path / "run_b"]
        for run in runs:
            run.mkdir()
            plot_n_per_sequence(scan, run / "n_count.png")
            plot_n_per_cycle_and_tile(scan, run / "n_cycle_tile.png")
            plot_individual_tiles(scan, run, 2)

        for name in [
            "n_count.png",
            "n_cycle_tile.png",
            "individual_tiles/x_R1_001_tile_1101.png",
            "individual_tiles/x_R1_001_tile_1102.png",
            "individual_tiles/x_R1_001_tile_1103.png",
        ]:
            assert (runs[0] / name).read_bytes() == (runs[1] / name).read_bytes(), name


# ---------------------------------------------------------------------------
# flags
# ---------------------------------------------------------------------------


class TestFlags:
    """Anomaly flags: bad tile, dark cycle and a single-cell outlier."""

    @pytest.fixture()
    def flagged(self, tmp_path):
        # 10 tiles x 40 reads, 20 cycles; tile 1110 has Ns everywhere,
        # cycle 5 is dark on every tile, tile 1103 has Ns at cycle 12 only
        records = []
        for t in range(1, 11):
            tile = f"11{t:02d}"
            for i in range(40):
                seq = list("ACGT" * 5)
                if tile == "1110" and i % 4 == 0:
                    seq[i % 20] = "N"
                if i < 30:
                    seq[4] = "N"
                if tile == "1103" and i < 20:
                    seq[11] = "N"
                records.append(_rec("".join(seq), tile=tile, x=str(i), y=str(i)))
        path = tmp_path / "S_S1_L001_R1_001.fastq.gz"
        write_fastq_gz(path, records)
        return build_tables(scan_files([path], "python"))

    def test_flag_kinds(self, flagged):
        flags = build_flags(flagged, warn_rate=0.01, fail_rate=0.05)
        got = {
            (r["kind"], r["tile"], r["cycle"], r["status"])
            for r in flags.iter_rows(named=True)
        }
        assert ("dark_cycle", None, 5, "FAIL") in got
        assert ("tile_cycle", 1103, 12, "FAIL") in got
        # 10 reads spread over 20 cycles: mean rate 0.0125 per cycle
        assert ("tile", 1110, None, "WARN") in got
        # the dark cycle is not repeated per tile and does not raise tile rates
        assert not any(k == "tile_cycle" and c == 5 for k, _, c, _ in got)
        # 1103's single bad cycle (20/40 reads) also lifts its mean to ~2.6%
        assert {t for k, t, _, _ in got if k == "tile"} == {1103, 1110}
        assert flags["status"].to_list() == sorted(
            flags["status"].to_list(), key=lambda s: s != "FAIL"
        )

    def test_high_thresholds_give_no_flags(self, run_files):
        tables = build_tables(scan_files(run_files, "python"))
        assert build_flags(tables, warn_rate=0.9, fail_rate=0.95).height == 0

    def test_sample_status(self, flagged):
        flags = build_flags(flagged)
        assert sample_status(flagged, flags).to_dicts() == [
            {"sample": "S", "status": "FAIL"}
        ]

    def test_report_and_multiqc_show_flags(self, flagged, tmp_path):
        flagged.flags = build_flags(flagged)
        scans = []  # summary tables only need the tables for flags
        page = write_report(scans, flagged, tmp_path).read_text()
        assert "Anomaly flags" in page and "dark_cycle" in page
        write_multiqc(flagged, tmp_path)
        table = json.loads(
            (tmp_path / "multiqc" / "nspector_flags_mqc.json").read_text()
        )
        assert (
            table["plot_type"] == "table" and len(table["data"]) == flagged.flags.height
        )
        stats = (tmp_path / "multiqc" / "nspector_general_stats_mqc.tsv").read_text()
        assert stats.splitlines()[-1].endswith("\tFAIL")


# ---------------------------------------------------------------------------
# InterOp tile metrics (--run-folder)
# ---------------------------------------------------------------------------


def _write_run_folder(run: pathlib.Path) -> None:
    """Write RunInfo.xml and TileMetricsOut.bin with the interop library itself."""
    np = pytest.importorskip("numpy")
    core = pytest.importorskip("interop.core")
    import interop.py_interop_metrics as metrics_mod
    import interop.py_interop_run as run_mod
    import interop.py_interop_run_metrics as run_metrics_mod

    (run / "InterOp").mkdir(parents=True)
    info = core._run_info_example_fixture(["1_1101", "1_1102"])
    metrics = run_metrics_mod.run_metrics(info)
    tiles = metrics.tile_metric_set()
    tiles.set_version(metrics_mod.tile_metric.LATEST_VERSION)
    for tile, pf in ((1101, 0.8), (1102, 0.5)):
        tiles.insert(metrics_mod.tile_metric(1, tile, 2e5, 2e5 * pf, 1e6, 1e6 * pf))
    info.write(str(run / "RunInfo.xml"))
    buf = np.zeros(metrics.calculate_buffer_size(run_mod.Tile), dtype=np.uint8)
    metrics.write_metrics_to_buffer(run_mod.Tile, buf)
    (run / "InterOp" / "TileMetricsOut.bin").write_bytes(buf.tobytes())


class TestRunMetrics:
    """--run-folder: InterOp tile metrics joined to N rates."""

    def test_read_tile_metrics(self, tmp_path):
        _write_run_folder(tmp_path / "run")
        frame = read_tile_metrics(tmp_path / "run")
        assert frame.select("read", "lane", "tile", "pct_pf").rows() == [
            ("R1", 1, 1101, 80.0),
            ("R1", 1, 1102, 50.0),
        ]
        # no Q metrics written: Q30 is null rather than 0
        assert frame["pct_q30"].null_count() == frame.height

    def test_missing_interop_package(self, tmp_path, monkeypatch):
        monkeypatch.setitem(sys.modules, "interop", None)
        with pytest.raises(RuntimeError, match=r"biomate\[interop\]"):
            read_tile_metrics(tmp_path)

    def test_join_and_report(self, run_files, tmp_path):
        tables = build_tables(scan_files(run_files, "python"))
        metrics = polars.DataFrame(
            [("R1", 1, 1101, 210.0, 80.0, 91.0)],
            schema=INTEROP_SCHEMA,
            orient="row",
        )
        tables.tile_interop = join_tile_interop(tables, metrics)
        row = tables.tile_interop.filter(
            (polars.col("lane") == 1) & (polars.col("tile") == 1101)
        )
        assert row["pct_pf"].item() == 80.0 and row["n_rate"].item() == 0.25
        # tiles without InterOp rows stay, with null metrics
        assert tables.tile_interop["pct_pf"].null_count() == 3
        spec = build_chart(tables).to_dict()
        assert "interop_metric" in {p["name"] for p in spec["params"]}
        pick = [p for p in spec["params"] if p["name"] == "pick_tile"]
        assert len(pick) == 1
        # one click handler each: tile map, tile x cycle heatmap, InterOp scatter
        vega = json.dumps(vl_convert.vegalite_to_vega(spec))
        assert vega.count('"name": "pick_tile_tuple", "on"') == 3
        assert "nspector_tile_interop.parquet" in {
            p.name for p in tables.write_parquet(tmp_path)
        }

    def test_run_folder_needs_interop_dir(self, tmp_path):
        fq = tmp_path / "x.fastq.gz"
        fq.touch()
        args = argparse.Namespace(
            input=[fq], output=tmp_path / "out", run_folder=tmp_path
        )
        with pytest.raises(argparse.ArgumentTypeError, match="InterOp"):
            validate_args(args)


def test_lowq_counter_batches_match_single_pass():
    """Flushing every few reads (many tiles, mixed lengths) gives the same counts."""
    rng = random.Random(5)
    reads = [
        (
            1,
            rng.choice((1101, 1102, 1103)),
            "".join(rng.choice("ACGTN") for _ in range(n)),
            "".join(chr(33 + rng.randint(0, 40)) for _ in range(n)),
        )
        for n in (rng.choice((4, 9)) for _ in range(500))
    ]
    small, big = LowQCounter(20, batch_size=3), LowQCounter(20, batch_size=10_000)
    for counter in (small, big):
        for read in reads:
            counter.add(*read)
    assert small.finish() == big.finish()
    assert small.reads == big.reads
    assert small.buffered == 0


# ---------------------------------------------------------------------------
# directory input and projects
# ---------------------------------------------------------------------------


def _demux_tree(root: pathlib.Path) -> dict[str, pathlib.Path]:
    """BCL Convert + bcl2fastq style tree with Undetermined at the root."""
    files = {
        "bclconvert": root / "ProjA" / "S1_S1_L001_R1_001.fastq.gz",
        "bcl2fastq": root / "ProjB" / "Sample_S2" / "S2_S2_L001_R1_001.fastq.gz",
        "undetermined": root / "Undetermined_S0_L001_R1_001.fastq.gz",
    }
    for i, path in enumerate(files.values()):
        path.parent.mkdir(parents=True, exist_ok=True)
        # all on tile 1101; 4 reads each, 1 + i of them with an N at cycle 2
        write_fastq_gz(
            path,
            [_rec("ANGT" if r <= i else "ACGT", x=str(r), y=str(r)) for r in range(4)],
        )
    (root / "ProjA" / "SampleSheet.csv").write_text("x")
    return files


class TestDirectoryInput:
    """Recursive discovery, project naming and project-aware tables."""

    def test_discover_and_project(self, tmp_path):
        files = _demux_tree(tmp_path)
        found = discover_fastqs([tmp_path, files["bclconvert"]])
        assert sorted(found) == sorted(files.values())  # duplicate dropped
        assert {k: project_of(p) for k, p in files.items()} == {
            "bclconvert": "ProjA",
            "bcl2fastq": "ProjB",
            "undetermined": "Undetermined",
        }

    def test_tables_split_by_project(self, tmp_path):
        _demux_tree(tmp_path)
        scans = scan_files(discover_fastqs([tmp_path]), "python")
        tables = build_tables(scans)
        assert set(tables.sample_lane["project"]) == {"ProjA", "ProjB", "Undetermined"}
        per_project = dict(
            tables.tile_cycle.filter(polars.col("cycle") == 2)
            .select("project", "n")
            .iter_rows()
        )
        assert per_project == {"ProjA": 1, "ProjB": 2, "Undetermined": 3}
        assert set(tables.tile_cycle["total"]) == {4}
        # flags stay flowcell-level: 6 of 12 reads with an N at cycle 2
        dense = dense_tile_cycle(tables).filter(polars.col("cycle") == 2)
        assert dense.select("n", "total").row(0) == (6, 12)
        anchored = with_tile_anchors(window_tile_cycle(tables)[0], tables)
        assert anchored.filter(polars.col("n") == 0).height == 3

    def test_rust_engine_matches_python_per_file_pool(self, tmp_path):
        pytest.importorskip("biomate_kernels")
        paths = discover_fastqs([_demux_tree(tmp_path)["bclconvert"].parents[1]])
        done = []
        rust = scan_files(paths, "rust", threads=2, done=done.append)
        assert rust == scan_files(paths, "python")
        assert sorted(done) == sorted(paths)


def test_read_control_offers_all_with_two_reads(tmp_path):
    paths = []
    for read in ("R1", "R2"):
        path = tmp_path / f"S_S1_L001_{read}_001.fastq.gz"
        write_fastq_gz(path, [_rec("ANGT"), _rec("ACGT", tile="1102")])
        paths.append(path)
    spec = build_chart(build_tables(scan_files(paths, "python"))).to_dict()
    read_sel = next(p for p in spec["params"] if p["name"] == "read_sel")
    assert read_sel["value"] == "All"
    assert read_sel["bind"]["options"] == ["All", "R1", "R2"]
    # legends toggle reads (shared by the read charts) and N / low-Q series
    legend = {p["name"]: p for p in spec["params"] if p.get("bind") == "legend"}
    assert set(legend) == {"pick_read", "pick_series"}
    assert len(legend["pick_read"]["views"]) >= 4


def test_project_summary_and_header(tmp_path):
    root = tmp_path / "run"
    _demux_tree(root)
    scans = scan_files(discover_fastqs([root]), "python")
    tables = build_tables(scans)
    tables.flags = build_flags(tables)
    summary = project_summary(tables)
    assert summary["project"].to_list() == ["ProjA", "ProjB", "Undetermined"]
    assert summary["reads with N"].to_list() == [1, 2, 3]
    assert set(summary["status"]) <= {"PASS", "WARN", "FAIL"}
    page = write_report(scans, tables, tmp_path).read_text()
    for text in (
        "<h2>Projects</h2>",
        "N per sequenced base",
        'class="badge',
        "window.print()",
        "prefers-color-scheme: dark",
        'id="theme-toggle"',
        "localStorage.setItem('nspector-theme'",
    ):
        assert text in page
    spec = json.dumps(build_chart(tables).to_dict())
    for section in ("Samples and cycles", "N within reads", "Flowcell"):
        assert section in spec


def test_robust_z_is_zero_at_the_median_of_a_zero_spread_group():
    """A tile at the median of a zero-MAD group scores 0, not NaN (>= any z)."""
    from biomate.nspector.flags import _robust_z

    df = polars.DataFrame({"g": ["a"] * 3, "rate": [0.02, 0.02, 0.0]})
    z = df.select(_robust_z(["g"]).alias("z")).get_column("z").to_list()
    assert z[:2] == [0.0, 0.0]
    assert z[2] == float("-inf")
