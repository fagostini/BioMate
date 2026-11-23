"""Tests for the barcode findings of triage."""

from biomate.triage.barcodes import classify, read_sample_counts, read_samples
from tests.conftest import DUAL_INDEX_SS


def _sample(sample_id, index1, index2=None, lane="L001", project="P"):
    return {
        "Lane": lane,
        "Sample_ID": sample_id,
        "Sample_Project": project,
        "index1": index1,
        "index2": index2,
    }


A = _sample("A", "AACCGGTA", "CCCAAATT")
B = _sample("B", "TTGACAGC", "GATCAGTA")


def _labels(samples, rows, counts=None):
    return [
        (f["Label"], f["Sample_ID"], f["Other_Sample_ID"])
        for f in classify(samples, rows, counts)
    ]


def _row(index, count=100, lane="L001"):
    return {"Lane": lane, "Undetermined_Index": index, "Count": count}


def test_index_mixing():
    assert _labels([A, B], [_row("AACCGGTA+GATCAGTA")]) == [("Index Mixing", "A", "B")]


def test_index_swap():
    assert _labels([A, B], [_row("CCCAAATT+AACCGGTA")]) == [("Index Swap", "A", None)]


def test_reverse_and_complement():
    assert _labels([A], [_row("ATGGCCAA+TTTTTTTT")]) == [("Reverse", "A", None)]
    assert _labels([A], [_row("TTGGCCAT+TTTTTTTT")]) == [("Complement", "A", None)]


def test_palindrome_and_own_barcode_are_not_findings():
    palindrome = _sample("C", "ACGTTGCA", "CCCAAATT")  # i7 reverse is itself
    assert _labels([palindrome], [_row("ACGTTGCA+TTTTTTTT")]) == []
    assert _labels([A, B], [_row("AACCGGTA+CCCAAATT")]) == []


def test_homopolymer_and_other_lane_dropped():
    assert _labels([A, B], [_row("AACCGGTA+GGGGAGTA")]) == []
    assert _labels([A, B], [_row("AACCGGTA+GATCAGTA", lane="L002")]) == []


def test_counts_threshold_and_putative_sample():
    counts = {("L001", "A"): 100, ("L001", "B"): 120}
    # Below the smallest sample: ignored
    assert _labels([A, B], [_row("AACCGGTA+GATCAGTA", 50)], counts) == []
    # Unexplained but sample-sized: putative sample of the lane's project
    found = classify([A, B], [_row("CATCATCA+TGTGTGTG", 110)], counts)
    assert [(f["Label"], f["Sample_Project"]) for f in found] == [
        ("Putative Sample", "P")
    ]
    # No counts: no putative sample
    assert _labels([A, B], [_row("CATCATCA+TGTGTGTG", 110)]) == []


def test_read_samples_and_counts(tmp_path):
    sheet = tmp_path / "SampleSheet.csv"
    sheet.write_text(DUAL_INDEX_SS)
    samples = read_samples(sheet)
    assert [(s["Lane"], s["Sample_ID"], s["index2"]) for s in samples] == [
        ("L001", "S1", "CGCTCCACGA"),
        ("L001", "S2", "TATCTTGTAG"),
        ("L002", "S3", "TTCCTATCAG"),
    ]
    assert read_sample_counts(tmp_path) is None
    reports = tmp_path / "Reports"
    reports.mkdir()
    (reports / "Demultiplex_Stats.csv").write_text(
        "Lane,SampleID,Sample_Project,Index,# Reads\n"
        "1,S1,ProjectA,GAACTGAGCG-CGCTCCACGA,10\n"
        "1,Undetermined,,,5\n"
    )
    assert read_sample_counts(tmp_path) == {
        ("L001", "S1"): 10,
        ("L001", "Undetermined"): 5,
    }
