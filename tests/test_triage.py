"""Tests for the triage module."""

import argparse
import json

import polars
import pytest

from biomate.triage.triage import main
from tests.conftest import DUAL_INDEX_SS, illumina_header, write_fastq_gz


def _record(lane: str, index: str) -> tuple[str, str, str]:
    return (illumina_header(lane=lane, index=index), "ACGT", "####")


@pytest.mark.parametrize("threads", [1, 2])
def test_main_writes_every_report(tmp_path, threads):
    """One run gives the strainer, top-index and waif reports."""
    run = tmp_path / "FC"
    run.mkdir()
    (run / "SampleSheet.csv").write_text(DUAL_INDEX_SS)
    # Lane 1: sample S1 with a reverse-complemented index2, plus a 1-mismatch variant
    write_fastq_gz(
        run / "Undetermined_S0_L001_R1_001.fastq.gz",
        [_record("1", "GAACTGAGCG+TCGTGGAGCG")] * 3
        + [_record("1", "GAACTGAGCT+TCGTGGAGCG")],
    )
    # Lane 2: sample S1 of lane 1 as written in the sheet
    write_fastq_gz(
        run / "Undetermined_S0_L002_R1_001.fastq.gz",
        [_record("2", "GAACTGAGCG+CGCTCCACGA")],
    )
    for dataset in ("RerunA", "RerunB"):
        stats = run / dataset / "Stats"
        stats.mkdir(parents=True)
        (stats / "Stats.json").write_text(
            json.dumps(
                {
                    "UnknownBarcodes": [
                        {"Lane": 1, "Barcodes": {"AAAACCCC+GGGGTTTT": 20}}
                    ]
                }
            )
        )
    out = tmp_path / f"out{threads}"
    out.mkdir()
    main(
        argparse.Namespace(
            input_path=run,
            output_path=out,
            sample_sheet=None,
            threads=threads,
            verbose=False,
            quiet=True,
        )
    )

    unexpected = polars.read_csv(out / "FC_unexpected_indexes.csv")
    assert sorted(
        unexpected.select("Lane", "Undetermined_Index", "Orientation").rows()
    ) == [
        ("L001", "GAACTGAGCG+TCGTGGAGCG", "index2_rc"),
        ("L001", "GAACTGAGCT+TCGTGGAGCG", "index2_rc"),
        ("L002", "GAACTGAGCG+CGCTCCACGA", "expected"),
    ]
    assert (out / "FC_lane_summary.csv").is_file()

    top = polars.read_csv(out / "FC_top_undetermined.csv")
    assert top.rows() == [
        ("Undetermined_S0_L001_R1_001.fastq.gz", "GAACTGAGCG+TCGTGGAGCG", 4, 1.0, 2),
        ("Undetermined_S0_L002_R1_001.fastq.gz", "GAACTGAGCG+CGCTCCACGA", 1, 1.0, 1),
    ]

    always = polars.read_csv(out / "FC_always_undetermined.csv")
    assert always.get_column("Barcode").to_list() == ["AAAACCCC+GGGGTTTT"]


def test_main_skips_waif_without_reruns(tmp_path):
    """Without Stats.json files only the FASTQ-based reports are written."""
    run = tmp_path / "FC"
    run.mkdir()
    (run / "SampleSheet.csv").write_text(DUAL_INDEX_SS)
    write_fastq_gz(
        run / "Undetermined_S0_L002_R1_001.fastq.gz",
        [_record("2", "TTTTTTTTTT+AAAAAAAAAA")],
    )
    main(
        argparse.Namespace(
            input_path=run,
            output_path=tmp_path,
            sample_sheet=None,
            threads=1,
            verbose=False,
            quiet=True,
        )
    )
    assert (tmp_path / "FC_top_undetermined.csv").is_file()
    assert not (tmp_path / "FC_always_undetermined.csv").exists()
    assert not (tmp_path / "FC_unexpected_indexes.csv").exists()
