"""Tests for the shared header index counter."""

import gzip

import pytest

from biomate.headers import count_header_indexes

RECORDS = "".join(
    f"@M:1:FC:1:1101:{i}:1 1:N:0:{idx}\nACGT\n+\n@III\n"
    for i, idx in enumerate(["ACGT+TTTT", "ACGT+TTTT", "GGCC+AATT"])
)


def test_plain_and_gzip_give_same_counts(tmp_path):
    """A quality line starting with '@' is not taken for a header."""
    plain, gz = tmp_path / "r.fastq", tmp_path / "r.fastq.gz"
    plain.write_text(RECORDS)
    with gzip.open(gz, "wt") as fh:
        fh.write(RECORDS)
    expected = {"ACGT+TTTT": 2, "GGCC+AATT": 1}
    assert count_header_indexes(plain) == count_header_indexes(gz) == expected


def test_empty_file_gives_no_counts(tmp_path):
    f = tmp_path / "r.fastq"
    f.write_text("")
    assert count_header_indexes(f) == {}


@pytest.mark.parametrize(
    "header", ["@M:1:FC:1:1101:1:1", "@read comment", ">r 1:N:0:A"]
)
def test_header_without_index_raises(tmp_path, header):
    f = tmp_path / "r.fastq"
    f.write_text(f"{header}\nACGT\n+\nIIII\n")
    with pytest.raises(ValueError, match="Invalid read name format"):
        count_header_indexes(f)
