"""Tests for the shared sample sheet reader."""

import pytest

from biomate.samplesheet import read_section


def test_read_section_tolerates_sheet_quirks(tmp_path):
    """Trailing commas, blank lines and a following section do not break parsing."""
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(
        "[Header],,,\nFileFormatVersion,2,,\n\n"
        "[BCLConvert_Data],,,\nLane,Sample_ID,index,\n1,S1,ACGT,\n\n2,S2,,\n"
        "[Cloud_Data],,,\nSample_ID,ProjectName,,\n"
    )
    data = read_section(sheet, infer_schema=False)
    assert data.columns == ["Lane", "Sample_ID", "index"]
    assert data.rows() == [("1", "S1", "ACGT"), ("2", "S2", None)]

    reads = tmp_path / "Reads.csv"
    reads.write_text("[Reads]\nRead1Cycles,151\nIndex1Cycles,10\n[Other]\n")
    assert read_section(reads, ("[Reads]",), has_header=False).rows() == [
        ("Read1Cycles", 151),
        ("Index1Cycles", 10),
    ]

    with pytest.raises(ValueError, match="Could not find section"):
        read_section(reads)
