"""Tests for the waif module."""

import argparse
import gzip
import json
import os
import re
import sys
from pathlib import Path

import polars as pl
import pytest

import biomate.waif.waif as waif
from biomate.waif.waif import (
    BoundedReader,
    _build_member_index,
    _count_for_file,
    _count_from_index,
    _count_undetermined,
    _find_member_offsets,
    _load_member_index,
    _member_index_for,
    get_shared_lanes_from_samplesheet,
    main,
    validate_args,
)

# Exact-mode tests fork the process pool while it may already be
# multi-threaded (polars workers started by earlier tests in this process).
# The children only run pure-Python gzip/Counter code and never touch
# polars, so the fork hazard flagged by CPython is theoretical here. In
# production runs the fork happens single-threaded (polars is first used
# after the pool work).
pytestmark = pytest.mark.filterwarnings(
    "ignore:This process .* multi-threaded, use of fork:DeprecationWarning"
)

LINUX_ONLY = pytest.mark.skipif(
    sys.platform != "linux", reason="requires os.fork (Linux only)"
)


def make_args(
    input_dir: Path,
    output_file: Path | None,
    exact: bool = False,
    sample_sheet: Path | None = None,
    block_cache: Path | None = None,
) -> argparse.Namespace:
    """Build the argparse.Namespace expected by validate_args and main."""
    return argparse.Namespace(
        input_path=input_dir,
        output_file=output_file,
        exact=exact,
        sample_sheet=sample_sheet,
        block_cache=block_cache,
        verbose=False,
        quiet=False,
    )


def run_main(
    input_dir: Path,
    output_file: Path,
    exact: bool = False,
    sample_sheet: Path | None = None,
    block_cache: Path | None = None,
) -> pl.DataFrame:
    """Run the waif module end to end and read the results CSV back."""
    main(make_args(input_dir, output_file, exact, sample_sheet, block_cache))
    return pl.read_csv(output_file)


def write_flowcell(root: Path, data: dict[str, dict[int, dict[str, int]]]) -> None:
    """Write ``data`` ({dataset: {lane: {barcode: count}}}) as Stats.json files."""
    for ds, lanes in data.items():
        p = root / ds / "Stats" / "Stats.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        unknown_barcodes = [
            {"Lane": lane, "Barcodes": barcodes} for lane, barcodes in lanes.items()
        ]
        p.write_text(json.dumps({"UnknownBarcodes": unknown_barcodes}))


def write_fastq(
    root: Path,
    dataset: str,
    lane: int,
    counts: dict[str, int],
    part: int = 1,
) -> None:
    """Write ``counts`` ({barcode: count}) as a gzipped R1 FASTQ for a lane."""
    f = root / dataset / f"Undetermined_S0_L{lane:03d}_R1_{part:03d}.fastq.gz"
    f.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(f, "wt") as out:
        for barcode, n in counts.items():
            for _ in range(n):
                out.write(f"@FC:1:1101:1:1101:1 1:N:0:{barcode}\nACGT\n+\nIIII\n")


def write_samplesheet(root: Path, lanes: set[int]) -> None:
    """Write a v2 SampleSheet where each lane carries two Sample_Projects."""
    lines = [
        "[Header]",
        "FileFormatVersion,2",
        "[BCLConvert_Data]",
        "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,"
        "OverrideCycles,BarcodeMismatchesIndex1,BarcodeMismatchesIndex2",
    ]
    for lane in sorted(lanes):
        lines.append(f"{lane},S_1,Sample_{lane}_1,AAAA,AAAA,ProjA,Y151;I4;I4;Y151,1,")
        lines.append(f"{lane},S_2,Sample_{lane}_2,CCCC,CCCC,ProjB,Y151;I4;I4;Y151,1,")
    (root / "SampleSheet.csv").write_text("\n".join(lines) + "\n")


class TestStatsJsonMode:
    """End-to-end tests of the default (Stats.json) mode."""

    def test_three_datasets_single_final_row(self, tmp_path: Path) -> None:
        """Only the full intersection is reported, once, as the clean k-way
        mean."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 200, "CCCC+CCCC": 200}},
                "D1": {1: {"AAAA+AAAA": 100, "WWWW+WWWW": 50}},
                "D2": {1: {"AAAA+AAAA": 300, "CCCC+CCCC": 300}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAA+AAAA"]
        assert out["Count"].to_list() == [200]  # round_up_10(mean(200, 100, 300))

    def test_mixed_lengths_prefix_aggregation(self, tmp_path: Path) -> None:
        """Two 10-base barcodes collapsing to one 8-base prefix have their
        counts summed before the mean; the result is in the shortest length
        space."""
        write_flowcell(
            tmp_path,
            {
                "D0": {
                    1: {
                        "AAAAAAAAAA+CCCCCCCCCC": 100,
                        "AAAAAAAATT+CCCCCCCCCC": 200,
                        "TTTTTTTTTT+AAAAAAAAAA": 400,
                    }
                },
                "D1": {
                    1: {
                        "AAAAAAAAAA+CCCCCCCCCC": 100,
                        "AAAAAAAATT+CCCCCCCCCC": 200,
                        "GGGGGGGGGG+TTTTTTTTTT": 50,
                    }
                },
                "D2": {1: {"AAAAAAAA+CCCCCCCC": 300, "TTTTTTTT+AAAAAAAA": 400}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAA+CCCCCCCC"]
        assert out["Count"].to_list() == [300]  # round_up_10(mean(300, 300, 300))

    def test_empty_intersection_no_crash(self, tmp_path: Path) -> None:
        """An empty running intersection must not crash the next step."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 10}},
                "D1": {1: {"CCCC+CCCC": 10}},
                "D2": {1: {"GGGG+GGGG": 10}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 0

    def test_k_way_mean_single_rounding(self, tmp_path: Path) -> None:
        """The k-way mean rounds once at the end (70), not per step (80)."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 1}},
                "D1": {1: {"AAAA+AAAA": 100}},
                "D2": {1: {"AAAA+AAAA": 100}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out["Count"].to_list() == [70]

    def test_mixed_length_pair_independent_of_names(self, tmp_path: Path) -> None:
        """A 10-base vs 8-base pair must match regardless of dataset names
        (the truncation direction comes from the length sort, not name/hash
        order)."""
        long_run = {1: {"AAAAAAAAAA+CCCCCCCCCC": 200}}
        short_run = {1: {"AAAAAAAA+CCCCCCCC": 100}}
        cases = [
            {"RunA": long_run, "RunZ": short_run},  # name order: long first
            {"RunA": short_run, "RunZ": long_run},  # name order: short first
        ]
        for i, data in enumerate(cases):
            root = tmp_path / f"case{i}"
            write_flowcell(root, data)
            out = run_main(root, root / "out.csv")
            assert out.height == 1, f"case {i}: expected a match"
            assert out["Barcode"].to_list() == ["AAAAAAAA+CCCCCCCC"]
            assert out["Count"].to_list() == [150]  # round_up_10(mean(200, 100))

    def test_dual_single_same_length(self, tmp_path: Path) -> None:
        """A dual-index and a single-index dataset on a lane compare on
        index1 only; the result is reported in the single-index form."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAAAAAAAA+CCCCCCCCCC": 200}},
                "D1": {1: {"AAAAAAAAAA": 100}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAAAA"]
        assert out["Count"].to_list() == [150]  # round_up_10(mean(200, 100))

    def test_dual_single_dual_i1_longer(self, tmp_path: Path) -> None:
        """When the dual dataset's index1 is longer, it is trimmed to the
        single-index length before comparing."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAAAAAAAAAA+CCCCCCCCCC": 200}},  # 12-base index1
                "D1": {1: {"AAAAAAAAAA": 100}},  # 10-base index1
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAAAA"]
        assert out["Count"].to_list() == [150]

    def test_dual_single_single_i1_longer(self, tmp_path: Path) -> None:
        """When the single-index dataset's index1 is longer, *it* is trimmed
        (the dual dataset is the hub by total length, so the other side must
        shrink)."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAAAAAA+CCCCCCCCCC": 200}},  # 8-base index1
                "D1": {1: {"AAAAAAAAAA": 100}},  # 10-base index1
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAA"]
        assert out["Count"].to_list() == [150]

    def test_single_single_mixed_length(self, tmp_path: Path) -> None:
        """Two single-index datasets of different lengths trim to the
        shorter."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAAAAAAAA": 200}},  # 10-base
                "D1": {1: {"AAAAAAAA": 100}},  # 8-base
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAA"]
        assert out["Count"].to_list() == [150]

    def test_mixed_dual_dual_single(self, tmp_path: Path) -> None:
        """A lane with two dual-index datasets and one single-index dataset is
        reported in the single-index form; the k-way mean uses all three
        counts."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAAAAAAAA+CCCCCCCCCC": 200}},
                "D1": {1: {"AAAAAAAAAA+CCCCCCCCCC": 100}},
                "D2": {1: {"AAAAAAAAAA": 300}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAAAAAAAA"]
        assert out["Count"].to_list() == [200]  # round_up_10(mean(200, 100, 300))

    def test_mixed_single_dual_within_one_dataset_is_an_error(
        self, tmp_path: Path
    ) -> None:
        """A single (Dataset, Lane) must be all-single or all-dual; mixing the
        two within one dataset is an error."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+BBBB": 10, "CCCC": 10}},
                "D1": {1: {"AAAA": 10}},
            },
        )
        with pytest.raises(ValueError, match="one i1/i2 index length"):
            run_main(tmp_path, tmp_path / "out.csv")

    def test_no_shared_lanes_empty_output(self, tmp_path: Path) -> None:
        """Lanes that are not shared between datasets yield an empty CSV."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 10}},
                "D1": {2: {"CCCC+CCCC": 10}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 0

    def test_no_stats_files_is_an_error(self, tmp_path: Path) -> None:
        """An input directory without any Stats.json files raises ValueError."""
        with pytest.raises(ValueError, match="Stats.json"):
            run_main(tmp_path, tmp_path / "out.csv")

    def test_barcode_with_empty_index1_is_an_error(self, tmp_path: Path) -> None:
        """A single-index barcode (no separator) is valid, but a barcode whose
        index1 is empty (e.g. '+AAAA') is still an error."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"+AAAA": 10}},
                "D1": {1: {"AAAA+AAAA": 10}},
            },
        )
        with pytest.raises(ValueError, match="empty index1"):
            run_main(tmp_path, tmp_path / "out.csv")

    def test_mixed_index_lengths_in_one_dataset_is_an_error(
        self, tmp_path: Path
    ) -> None:
        """Mixed index lengths within one (Dataset, Lane) raise ValueError."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 10, "GGGGGGGGGG+AAAA": 10}},
                "D1": {1: {"AAAA+AAAA": 10}},
            },
        )
        with pytest.raises(ValueError, match="one i1/i2 index length"):
            run_main(tmp_path, tmp_path / "out.csv")

    def test_missing_unknown_barcodes_skipped(self, tmp_path: Path) -> None:
        """A Stats.json without UnknownBarcodes is skipped with a warning."""
        write_flowcell(tmp_path, {"D0": {1: {"AAAA+AAAA": 10}}})
        stray = tmp_path / "D1" / "Stats"
        stray.mkdir(parents=True)
        (stray / "Stats.json").write_text(json.dumps({"Flowcell": "X"}))
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 0  # only one usable dataset -> no shared lane

    def test_demultiplexing_dir_excluded(self, tmp_path: Path) -> None:
        """Files under Demultiplexing/ are ignored; only re-runs are
        compared."""
        write_flowcell(
            tmp_path,
            {
                "Demultiplexing": {1: {"AAAA+AAAA": 999999}},
                "Demultiplexing_0": {1: {"AAAA+AAAA": 200}},
                "Demultiplexing_1": {1: {"AAAA+AAAA": 200}},
            },
        )
        out = run_main(tmp_path, tmp_path / "out.csv")
        assert out.height == 1
        assert out["Count"].to_list() == [200]  # mean of the two re-runs only


class TestExactMode:
    """End-to-end tests of the --exact (FASTQ counting) mode."""

    def test_exact_uses_all_datasets(self, tmp_path: Path) -> None:
        """Every dataset's undetermined FASTQ must feed the intersection (a
        dict keyed by lane silently dropped all but the last file per
        lane)."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 200, "CCCC+CCCC": 100})
        write_fastq(tmp_path, "D1", 1, {"AAAA+AAAA": 100, "GGGG+GGGG": 50})
        write_samplesheet(tmp_path, {1})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAA+AAAA"]
        assert out["Count"].to_list() == [150]  # round_up_10(mean(200, 100))

    def test_exact_single_index(self, tmp_path: Path) -> None:
        """Exact mode counts bare single-index barcodes from the FASTQ headers
        and intersects them with a dual-index dataset on index1 only (single
        form)."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+BBBB": 200})
        write_fastq(tmp_path, "D1", 1, {"AAAA": 100})
        write_samplesheet(tmp_path, {1})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        assert out.height == 1
        assert out["Barcode"].to_list() == ["AAAA"]
        assert out["Count"].to_list() == [150]  # round_up_10(mean(200, 100))

    def test_exact_demultiplexing_dir_excluded(self, tmp_path: Path) -> None:
        """FASTQs under Demultiplexing/ are primary demux output, not a
        re-run."""
        write_fastq(tmp_path, "Demultiplexing", 1, {"AAAA+AAAA": 999999})
        write_fastq(tmp_path, "Demultiplexing_0", 1, {"AAAA+AAAA": 200})
        write_fastq(tmp_path, "Demultiplexing_1", 1, {"AAAA+AAAA": 200})
        write_samplesheet(tmp_path, {1})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        assert out.height == 1
        assert out["Count"].to_list() == [200]  # mean of the two re-runs only

    def test_exact_multi_part_files_merged(self, tmp_path: Path) -> None:
        """bcl-convert's 1M-read parts (_001, _002) of one dataset must be
        summed, not joined as separate datasets."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 100}, part=1)
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 100, "CCCC+CCCC": 50}, part=2)
        write_fastq(tmp_path, "D1", 1, {"AAAA+AAAA": 300, "CCCC+CCCC": 50})
        write_samplesheet(tmp_path, {1})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        assert out.height == 2
        results = dict(zip(out["Barcode"].to_list(), out["Count"].to_list()))
        assert results["AAAA+AAAA"] == 250  # round_up_10(mean(200, 300))
        assert results["CCCC+CCCC"] == 50  # round_up_10(mean(50, 50))

    def test_exact_lane_without_fastq_is_skipped(self, tmp_path: Path) -> None:
        """A samplesheet-selected lane with no FASTQ files must not crash the
        run."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 100})
        write_fastq(tmp_path, "D1", 1, {"AAAA+AAAA": 100})
        write_samplesheet(tmp_path, {1, 2})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        assert out.height == 1
        assert out["Lane"].to_list() == [1]
        assert out["Count"].to_list() == [100]

    def test_exact_single_file_serial_path(self, tmp_path: Path) -> None:
        """A single FASTQ file takes the serial path (no pool) and reports the
        lane's own counts, each rounded up to a multiple of 10."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 105, "CCCC+CCCC": 1})
        write_samplesheet(tmp_path, {1})
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True)
        results = dict(zip(out["Barcode"].to_list(), out["Count"].to_list()))
        assert results == {"AAAA+AAAA": 110, "CCCC+CCCC": 10}

    def test_exact_output_is_deterministic(self, tmp_path: Path) -> None:
        """The parallel extraction must not depend on file processing
        order."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 100, "CCCC+CCCC": 50})
        write_fastq(
            tmp_path, "D1", 1, {"AAAA+AAAA": 300, "CCCC+CCCC": 50, "GGGG+GGGG": 10}
        )
        write_samplesheet(tmp_path, {1})
        run_main(tmp_path, tmp_path / "out1.csv", exact=True)
        run_main(tmp_path, tmp_path / "out2.csv", exact=True)
        assert (tmp_path / "out1.csv").read_bytes() == (
            tmp_path / "out2.csv"
        ).read_bytes()


def _reference_counts(text: str) -> dict[str, int]:
    """Count the index pairs of the FASTQ header lines, line by line."""
    counts: dict[str, int] = {}
    for line in text.split("\n"):
        if line.startswith("@"):
            idx = line.rsplit(":", 1)[-1].rstrip()
            counts[idx] = counts.get(idx, 0) + 1
    return counts


def _to_bases(n: int, width: int) -> str:
    """Encode the bits of n as a DNA sequence of the given width."""
    return "".join("ACGT"[(n >> (2 * j)) & 3] for j in range(width))


def fastq_fixtures(tmp_path: Path) -> dict[str, tuple[Path, dict[str, int]]]:
    """Write gzipped FASTQ variants and return {name: (path, expected counts)}.

    Covers CRLF endings, a missing final newline, and non-header lines
    containing '@' and ':' characters.
    """
    lines = []
    for i in range(200):
        barcode = f"{_to_bases(i, 10)}+{_to_bases(3 * i + 1, 10)}"
        lines.append(f"@FC:1:1101:1:{1000 + i}:1 1:N:0:{barcode}")
        lines.append("ACGTACGT")
        lines.append("+")
        lines.append("I@I:IIII:ACGT+I#")  # quality line with '@' and ':'
    text = "\n".join(lines) + "\n"

    fixtures = {
        "lf": text,
        "crlf": "\r\n".join(lines) + "\r\n",
        # Final header without a trailing newline (malformed file).
        "no-final-nl": "\n".join(lines + ["@FC:1:1101:1:9999:9 1:N:0:ZZZZ+QQQQ"]),
    }
    out: dict[str, tuple[Path, dict[str, int]]] = {}
    for name, content in fixtures.items():
        f = tmp_path / f"test_{name}.fastq.gz"
        with gzip.open(f, "wb") as fh:
            fh.write(content.encode())
        out[name] = (f, _reference_counts(content))
    return out


class TestCounters:
    """The counting passes must equal a line-by-line reference."""

    def test_counter_matches_line_by_line_reference(self, tmp_path: Path) -> None:
        """The serial counter must equal the line-by-line reference."""
        for name, (f, expected) in fastq_fixtures(tmp_path).items():
            got = dict(_count_undetermined(f))
            assert got == expected, f"{name}"

    def test_threaded_counter_matches_line_by_line_reference(
        self, tmp_path: Path
    ) -> None:
        """The threaded inflate+count pipeline must equal the same
        reference."""
        for name, (f, expected) in fastq_fixtures(tmp_path).items():
            got = dict(_count_undetermined(f, force_threaded=True))
            assert got == expected, f"{name}"

    def test_threaded_multi_member_mid_chunk(self, tmp_path: Path) -> None:
        """A multi-member source (boundary mid-chunk) must be counted
        fully."""
        text = _fastq_text(64, 0) + _fastq_text(64, 100)
        f = tmp_path / "midchunk.fastq.gz"
        f.write_bytes(
            gzip.compress(_fastq_text(64, 0).encode())
            + gzip.compress(_fastq_text(64, 100).encode())
        )
        got = dict(_count_undetermined(f, force_threaded=True))
        assert got == _reference_counts(text)

    def test_threaded_multi_member_three(self, tmp_path: Path) -> None:
        """Three concatenated members must all be counted."""
        parts = [_fastq_text(40, 10 * k) for k in range(3)]
        f = tmp_path / "three.fastq.gz"
        f.write_bytes(b"".join(gzip.compress(p.encode()) for p in parts))
        got = dict(_count_undetermined(f, force_threaded=True))
        assert got == _reference_counts("".join(parts))

    def test_threaded_multi_member_aligned_boundary(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """Regression: a member ending exactly at a chunk edge must not cause
        later members to be dropped (the producer used to stop at the
        edge)."""
        m1 = gzip.compress(_fastq_text(50, 0).encode())
        m2 = gzip.compress(_fastq_text(50, 100).encode())
        f = tmp_path / "aligned.fastq.gz"
        f.write_bytes(m1 + m2)
        monkeypatch.setattr(waif, "_INFLATE_CHUNK", len(m1))
        got = dict(_count_undetermined(f, force_threaded=True))
        assert got == _reference_counts(_fastq_text(50, 0) + _fastq_text(50, 100))

    def test_threaded_member_spans_multiple_chunks(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A member larger than the chunk (spanning many reads) followed by a
        second member must both be counted fully (the production-scale
        case)."""
        p1 = _fastq_text(500, 0)
        p2 = _fastq_text(100, 1000)
        f = tmp_path / "bigmember.fastq.gz"
        f.write_bytes(gzip.compress(p1.encode()) + gzip.compress(p2.encode()))
        monkeypatch.setattr(waif, "_INFLATE_CHUNK", 64)
        got = dict(_count_undetermined(f, force_threaded=True))
        assert got == _reference_counts(p1 + p2)


def _fastq_text(n: int, offset: int) -> str:
    """Return ``n`` FASTQ records with barcodes derived from ``i + offset``."""
    lines: list[str] = []
    for i in range(n):
        j = i + offset
        barcode = f"{_to_bases(j, 10)}+{_to_bases(3 * j + 1, 10)}"
        lines.append(f"@FC:1:1101:1:{j}:1 1:N:0:{barcode}")
        lines.extend(["ACGT", "+", "IIII"])
    return "\n".join(lines) + "\n"


def _aligned_multi_member(
    tmp_path: Path, name: str, texts: list[str]
) -> tuple[Path, dict[str, int]]:
    """A multi-member gzip file whose member boundaries fall on line
    boundaries (each member is a whole number of FASTQ records)."""
    f = tmp_path / name
    with open(f, "wb") as out:
        for t in texts:
            out.write(gzip.compress(t.encode()))
    return f, _reference_counts("".join(texts))


def _misaligned_two_member(
    tmp_path: Path, name: str, text: str
) -> tuple[Path, dict[str, int]]:
    """A two-member gzip file whose single boundary falls inside a
    header."""
    enc = text.encode()
    starts = [m.start() for m in re.finditer(rb"@", enc)]
    cut = starts[len(starts) // 2] + 5  # mid-header, never on a line boundary
    f = tmp_path / name
    with open(f, "wb") as out:
        out.write(gzip.compress(enc[:cut]))
        out.write(gzip.compress(enc[cut:]))
    return f, _reference_counts(text)


class TestMemberIndex:
    """Gzip member-offset index and the forked parallel count."""

    def test_member_index_build_and_count(self, tmp_path: Path) -> None:
        """Building the index and counting must equal the serial reference,
        with valid metadata and a line-aligned flag."""
        for name, (f, expected) in fastq_fixtures(tmp_path).items():
            cache_dir = tmp_path / "idx"
            index_path = _member_index_for(f, cache_dir)
            index = _build_member_index(f, index_path)
            assert index["src"] == str(f.resolve())
            assert index["offsets"][0] == 0
            assert index["offsets"][-1] == f.stat().st_size
            assert index["mtime_ns"] == f.stat().st_mtime_ns
            assert index["size"] == f.stat().st_size
            assert index["aligned"] is True
            # Single member: the bounded range covers the whole file.
            got = dict(_count_from_index(f, index, budget=3))
            assert got == expected, f"{name}"

    @LINUX_ONLY
    def test_count_from_index_parallel_aligned(self, tmp_path: Path) -> None:
        """A line-aligned multi-member file counted by forked workers over
        member ranges must equal the serial reference."""
        f, expected = _aligned_multi_member(
            tmp_path, "aligned.fastq.gz", [_fastq_text(40, 10 * k) for k in range(4)]
        )
        index = _build_member_index(f, _member_index_for(f, tmp_path / "idx"))
        assert index["aligned"] is True
        assert len(index["offsets"]) - 1 == 4
        got = dict(_count_from_index(f, index, budget=4))
        assert got == expected

    @LINUX_ONLY
    def test_count_for_file_uses_index_when_aligned(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """_count_for_file must route a large line-aligned multi-member file
        to the index (forked) path; a second call hits the stored
        index."""
        monkeypatch.setattr(waif, "THREADED_MIN_SIZE", 0)
        f, expected = _aligned_multi_member(
            tmp_path, "aligned.fastq.gz", [_fastq_text(40, 10 * k) for k in range(3)]
        )
        cache_dir = tmp_path / "idx"
        calls: list[int] = []
        real = _count_from_index

        def spy(ff: Path, index: dict, budget: int):
            """Record that the index path was taken, then delegate to the real counter."""
            calls.append(1)
            return real(ff, index, budget)

        monkeypatch.setattr(waif, "_count_from_index", spy)
        got = dict(_count_for_file(f, budget=3, cache_dir=cache_dir))
        assert got == expected
        assert calls, "expected the index (parallel) path to be taken"
        assert _load_member_index(f, _member_index_for(f, cache_dir)) is not None
        got2 = dict(_count_for_file(f, budget=3, cache_dir=cache_dir))
        assert got2 == expected

    @LINUX_ONLY
    def test_count_for_file_falls_back_when_misaligned(
        self, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """A file whose member boundary falls mid-header must NOT use the
        index path; it must fall back to the (correct) serial count."""
        monkeypatch.setattr(waif, "THREADED_MIN_SIZE", 0)
        f, expected = _misaligned_two_member(
            tmp_path, "misaligned.fastq.gz", _fastq_text(60, 0)
        )
        _, aligned = _find_member_offsets(f)
        assert aligned is False
        index = _build_member_index(f, _member_index_for(f, tmp_path / "idx"))
        assert index["aligned"] is False
        # Document that the unguarded index path would be wrong here.
        assert dict(_count_from_index(f, index, budget=2)) != expected
        got = dict(_count_for_file(f, budget=3, cache_dir=tmp_path / "idx"))
        assert got == expected  # correct via the serial fallback

    def test_index_invalidation(self, tmp_path: Path) -> None:
        """A source mtime change invalidates the index; a rebuild repairs
        it."""
        f, expected = fastq_fixtures(tmp_path)["lf"]
        cache_dir = tmp_path / "idx"
        index_path = _member_index_for(f, cache_dir)
        _build_member_index(f, index_path)
        assert _load_member_index(f, index_path) is not None
        st = f.stat()
        os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns + 1))
        assert _load_member_index(f, index_path) is None
        _build_member_index(f, index_path)
        index = _load_member_index(f, index_path)
        assert index is not None
        got = dict(_count_from_index(f, index, budget=2))
        assert got == expected

    def test_bounded_reader_slice(self, tmp_path: Path) -> None:
        """BoundedReader over a member range must expose exactly those
        members' decompressed bytes and stop at the range end."""
        parts = [_fastq_text(20, 100 * k) for k in range(3)]
        f = tmp_path / "mm.fastq.gz"
        with open(f, "wb") as out:
            for p in parts:
                out.write(gzip.compress(p.encode()))
        index = _build_member_index(f, _member_index_for(f, tmp_path / "idx"))
        offsets = index["offsets"]
        with gzip.GzipFile(
            fileobj=BoundedReader(f, offsets[1], offsets[2]), mode="rb"
        ) as fh:
            assert fh.read().decode() == parts[1]


class TestSampleSheet:
    """Tests for get_shared_lanes_from_samplesheet."""

    def test_shared_lanes_from_v1_samplesheet(self, tmp_path: Path) -> None:
        """FCID-format sheets read Sample_Project from the 12th field; rows
        with an empty project are ignored."""
        sheet = tmp_path / "SampleSheet.csv"
        sheet.write_text(
            "FCID,Lane,Sample_ID,Sample_Name,Sample_Ref,index,index2,Description,"
            "Control,Recipe,Operator,Sample_Project\n"
            "X,1,S_1,N_1,Ref,AAAA,AAAA,d,N,28-90,Op,ProjA\n"
            "X,1,S_2,N_2,Ref,CCCC,CCCC,d,N,28-90,Op,ProjB\n"
            "X,2,S_3,N_3,Ref,GGGG,GGGG,d,N,28-90,Op,ProjC\n"
            "X,2,S_4,N_4,Ref,NOINDEX,,d,N,28-90,Op,ProjC\n"
            "X,2,S_5,N_5,Ref,TTTT,TTTT,d,N,28-90,Op,\n"
        )
        assert get_shared_lanes_from_samplesheet(sheet) == [1]


class TestValidateArgs:
    """Tests for validate_args."""

    def _make_args(self, tmp_path: Path, **overrides) -> argparse.Namespace:
        """Build a Namespace with sensible defaults, overridable per test."""
        base = dict(
            input_path=tmp_path / "input",
            output_file=None,
            exact=False,
            sample_sheet=None,
            block_cache=None,
            verbose=False,
            quiet=False,
        )
        base.update(overrides)
        return argparse.Namespace(**base)

    def test_valid_args_returned(self, tmp_path: Path) -> None:
        """validate_args returns the same Namespace when all arguments are valid."""
        (tmp_path / "input").mkdir()
        args = self._make_args(tmp_path)
        assert validate_args(args) is args

    def test_invalid_input_path_raises(self, tmp_path: Path) -> None:
        """A non-existent input path raises ArgumentTypeError."""
        args = self._make_args(tmp_path)
        with pytest.raises(argparse.ArgumentTypeError, match="Input path"):
            validate_args(args)

    def test_output_dir_created_if_missing(self, tmp_path: Path) -> None:
        """The parent directory of the output file is created automatically."""
        (tmp_path / "input").mkdir()
        args = self._make_args(tmp_path, output_file=tmp_path / "new" / "out.csv")
        validate_args(args)
        assert (tmp_path / "new").is_dir()

    def test_exact_requires_sample_sheet(self, tmp_path: Path) -> None:
        """--exact without any Sample Sheet raises ArgumentTypeError."""
        (tmp_path / "input").mkdir()
        args = self._make_args(tmp_path, exact=True)
        with pytest.raises(argparse.ArgumentTypeError, match="Sample Sheet"):
            validate_args(args)

    def test_explicit_sample_sheet_validated(self, tmp_path: Path) -> None:
        """A non-existent explicit Sample Sheet raises ArgumentTypeError."""
        (tmp_path / "input").mkdir()
        args = self._make_args(
            tmp_path, exact=True, sample_sheet=tmp_path / "nonexistent.csv"
        )
        with pytest.raises(argparse.ArgumentTypeError, match="does not exist"):
            validate_args(args)

    def test_explicit_sample_sheet_accepted(self, tmp_path: Path) -> None:
        """An existing explicit Sample Sheet is accepted."""
        (tmp_path / "input").mkdir()
        sheet = tmp_path / "sheet.csv"
        sheet.write_text(
            "[BCLConvert_Data]\n"
            "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,"
            "OverrideCycles\n"
            "1,S_1,Sample_1_1,AAAA,AAAA,ProjA,Y151;I4;I4;Y151\n"
        )
        args = self._make_args(tmp_path, exact=True, sample_sheet=sheet)
        assert validate_args(args) is args

    @LINUX_ONLY
    def test_block_cache_dir_created(self, tmp_path: Path) -> None:
        """The --block-cache directory is created automatically."""
        (tmp_path / "input").mkdir()
        args = self._make_args(tmp_path, block_cache=tmp_path / "cache")
        validate_args(args)
        assert (tmp_path / "cache").is_dir()

    @pytest.mark.skipif(sys.platform == "linux", reason="guard only fires on non-Linux")
    def test_block_cache_rejected_on_non_linux(self, tmp_path: Path) -> None:
        """--block-cache is rejected on non-Linux platforms."""
        (tmp_path / "input").mkdir()
        args = self._make_args(tmp_path, block_cache=tmp_path / "cache")
        with pytest.raises(argparse.ArgumentTypeError, match="Linux"):
            validate_args(args)


class TestMainIntegration:
    """Behaviour of the main() entry point (sample-sheet option, stdout)."""

    def test_explicit_sample_sheet_option(self, tmp_path: Path) -> None:
        """--sample-sheet may point outside the input directory."""
        write_fastq(tmp_path, "D0", 1, {"AAAA+AAAA": 200})
        write_fastq(tmp_path, "D1", 1, {"AAAA+AAAA": 100})
        sheet = tmp_path / "elsewhere" / "Sheet.csv"
        sheet.parent.mkdir(parents=True)
        sheet.write_text(
            "[BCLConvert_Data]\n"
            "Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,"
            "OverrideCycles\n"
            "1,S_1,Sample_1_1,AAAA,AAAA,ProjA,Y151;I4;I4;Y151\n"
            "1,S_2,Sample_1_2,CCCC,CCCC,ProjB,Y151;I4;I4;Y151\n"
        )
        out = run_main(tmp_path, tmp_path / "out.csv", exact=True, sample_sheet=sheet)
        assert out["Barcode"].to_list() == ["AAAA+AAAA"]
        assert out["Count"].to_list() == [150]

    def test_results_written_to_stdout(self, tmp_path: Path, capsys) -> None:
        """Without --output-file the CSV goes to stdout."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 200}},
                "D1": {1: {"AAAA+AAAA": 100}},
            },
        )
        main(make_args(tmp_path, None))
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert lines[0] == "Lane,Barcode,Count,Total %"
        assert lines[1] == "1,AAAA+AAAA,150,100.0"

    def test_stdout_empty_csv_when_no_shared_lanes(
        self, tmp_path: Path, capsys
    ) -> None:
        """Without shared lanes the stdout CSV carries only the header."""
        write_flowcell(
            tmp_path,
            {
                "D0": {1: {"AAAA+AAAA": 10}},
                "D1": {2: {"CCCC+CCCC": 10}},
            },
        )
        main(make_args(tmp_path, None))
        captured = capsys.readouterr()
        assert captured.out.strip() == "Lane,Barcode,Count,Total %"
