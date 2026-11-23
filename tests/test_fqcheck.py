"""Tests for the fqcheck module."""

import argparse
import hashlib

import polars
import pytest

from biomate.fqcheck.fqcheck import main
from tests.conftest import illumina_header, write_fastq_gz


def _records(n: int, read: str = "1") -> list[tuple[str, str, str]]:
    return [
        (illumina_header(x=str(i)).replace(" 1:", f" {read}:"), "ACGT", "####")
        for i in range(n)
    ]


def _run(input_path, output_path, md5=False):
    main(
        argparse.Namespace(
            input_path=input_path,
            output_path=output_path,
            md5=md5,
            threads=1,
            verbose=False,
            quiet=True,
        )
    )
    return polars.read_csv(output_path / f"{input_path.name}_fqcheck.csv")


@pytest.fixture
def run(tmp_path):
    run = tmp_path / "FC"
    run.mkdir()
    for read in ("R1", "R2", "I1"):
        write_fastq_gz(run / f"A_S1_L001_{read}_001.fastq.gz", _records(3, read[1]))
    return run


def test_valid_set_is_ok_and_md5_matches(run, tmp_path):
    report = _run(run, tmp_path, md5=True)
    assert report.rows() == [
        (
            "A_S1_L001_*_001.fastq.gz",
            "A_S1_L001_R1_001.fastq.gz;A_S1_L001_R2_001.fastq.gz;A_S1_L001_I1_001.fastq.gz",
            3,
            "ok",
            None,
        )
    ]
    lines = (tmp_path / "FC.md5").read_text().splitlines()
    assert len(lines) == 3
    for line in lines:
        digest, name = line.split("  ")
        assert digest == hashlib.md5((run / name).read_bytes()).hexdigest()


@pytest.mark.parametrize(
    "fault",
    ["truncated", "short_mate", "name_out_of_sync", "short_quality", "no_r1"],
)
def test_faults_are_errors(run, tmp_path, fault):
    r1 = run / "A_S1_L001_R1_001.fastq.gz"
    r2 = run / "A_S1_L001_R2_001.fastq.gz"
    if fault == "truncated":
        r2.write_bytes(r2.read_bytes()[:-20])
    elif fault == "short_mate":
        write_fastq_gz(r2, _records(2, "2"))
    elif fault == "name_out_of_sync":
        records = _records(3, "2")
        records[1] = ("OTHER:1:AABCCC:1:1101:1:1 2:N:0:A", "ACGT", "####")
        write_fastq_gz(r2, records)
    elif fault == "short_quality":
        records = _records(3, "2")
        records[2] = (records[2][0], "ACGT", "###")
        write_fastq_gz(r2, records)
    elif fault == "no_r1":
        r1.unlink()
    with pytest.raises(SystemExit) as exit_info:
        _run(run, tmp_path)
    assert exit_info.value.code == 1
    report = polars.read_csv(tmp_path / "FC_fqcheck.csv")
    assert report["Status"].to_list() == ["error"]


def test_non_illumina_name_is_warning(tmp_path):
    run = tmp_path / "FC"
    run.mkdir()
    write_fastq_gz(run / "reads.fastq.gz", _records(2))
    report = _run(run, tmp_path)
    assert report["Status"].to_list() == ["warning"]
    assert "not an Illumina file name" in report["Message"][0]
