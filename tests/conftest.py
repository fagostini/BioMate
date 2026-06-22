"""Shared fixtures and helpers for BioMate tests."""

import gzip
import pathlib

import pytest

# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------


def write_fastq_gz(path: pathlib.Path, records: list[tuple[str, str, str]]) -> None:
    """Write a list of (header, seq, qual) tuples to a gzip-compressed FASTQ."""
    with gzip.open(path, "wt") as fh:
        for header, seq, qual in records:
            fh.write(f"@{header}\n{seq}\n+\n{qual}\n")


def illumina_header(
    tile: str = "1101",
    x: str = "10000",
    y: str = "20000",
    lane: str = "1",
    index: str = "ACGTACGT+TGCATGCA",
) -> str:
    """Return a minimal Illumina FASTQ read header string (without leading @)."""
    return f"VH00001:1:AABCCC:{lane}:{tile}:{x}:{y} 1:N:0:{index}"


# ---------------------------------------------------------------------------
# Sample sheet content strings
# ---------------------------------------------------------------------------

DUAL_INDEX_SS = """\
[Header]
FileFormatVersion,2
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
1,S1,Sample1,GAACTGAGCG,CGCTCCACGA,ProjectA,Y151;I10;I10;Y151
1,S2,Sample2,AGGTCAGATA,TATCTTGTAG,ProjectA,Y151;I10;I10;Y151
2,S3,Sample3,TAAACCCTAG,TTCCTATCAG,ProjectB,Y151;I10;I10;Y151
"""

SINGLE_INDEX_SS = """\
[Header]
FileFormatVersion,2
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,Sample_Project,OverrideCycles
1,S1,Sample1,GAACTGAGCG,ProjectA,Y151;I10;Y151
2,S2,Sample2,AGGTCAGATA,ProjectB,Y151;I10;Y151
"""

NO_DATA_SECTION_SS = """\
[Header]
FileFormatVersion,2
Lane,Sample_ID,index
1,S1,ACGT
"""

MISSING_COLUMNS_SS = """\
[Header]
FileFormatVersion,2
[BCLConvert_Data]
Sample_ID,Sample_Name,index
S1,Sample1,GAACTGAGCG
"""

EMPTY_DATA_SS = """\
[Header]
FileFormatVersion,2
[BCLConvert_Data]
Lane,Sample_ID,Sample_Name,index,index2,Sample_Project,OverrideCycles
"""


# ---------------------------------------------------------------------------
# pytest fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def dual_index_sample_sheet(tmp_path: pathlib.Path) -> pathlib.Path:
    path = tmp_path / "SampleSheet_dual.csv"
    path.write_text(DUAL_INDEX_SS)
    return path


@pytest.fixture
def single_index_sample_sheet(tmp_path: pathlib.Path) -> pathlib.Path:
    path = tmp_path / "SampleSheet_single.csv"
    path.write_text(SINGLE_INDEX_SS)
    return path


@pytest.fixture
def undetermined_fastq_l001(tmp_path: pathlib.Path) -> pathlib.Path:
    """A valid gzip FASTQ undetermined file for lane L001."""
    path = tmp_path / "Undetermined_S0_L001_R1_001.fastq.gz"
    records = [
        (illumina_header(lane="1", index="ACGTACGT+TGCATGCA"), "ACGT", "####"),
        (illumina_header(lane="1", index="ACGTACGT+TGCATGCA"), "TTTT", "$$$$"),
        (illumina_header(lane="1", index="CCCCCCCC+GGGGGGGG"), "AAAA", "@@@@"),
    ]
    write_fastq_gz(path, records)
    return path
