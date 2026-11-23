"""Explain undetermined barcodes with the sample sheet indexes of their lane.

Ported from flowceller's ``check_barcodes``, without the reverse-complement
rule (strainer reports it): Index Mixing, Index Swap, Reverse, Complement
and Putative Sample.
"""

import csv
import logging
import pathlib
import statistics

from biomate.samplesheet import read_section
from biomate.strainer.strainer import _COMPLEMENT

# A run of four G or N in either index is a no-signal artefact
DISCARD_RUNS = ("GGGG", "NNNN")
# Putative samples must be within this many stdevs of the project median
PUTATIVE_STDEVS = 1.5
FIELDS = [
    "Lane",
    "Label",
    "Undetermined_Index",
    "Count",
    "Sample_ID",
    "Other_Sample_ID",
    "Sample_Project",
    "Matched",
]


def _lane(value: str) -> str:
    """Lane token of the undetermined rows, e.g. ``1`` -> ``L001``."""
    return f"L{int(value):03d}"


def _rc(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def read_samples(sample_sheet: pathlib.Path) -> list[dict]:
    """Sheet rows with an index: Lane, Sample_ID, Sample_Project, index1, index2."""
    data = read_section(sample_sheet, infer_schema=False)
    missing = {"Lane", "Sample_ID", "index"} - set(data.columns)
    if missing:
        raise ValueError(f"Sample Sheet missing columns: {', '.join(sorted(missing))}")
    samples = []
    for row in data.iter_rows(named=True):
        index1 = (row["index"] or "").strip().upper()
        if not index1:
            continue
        index2 = (row.get("index2") or "").strip().upper() or None
        samples.append(
            {
                "Lane": _lane(row["Lane"]),
                "Sample_ID": row["Sample_ID"],
                "Sample_Project": row.get("Sample_Project"),
                "index1": index1,
                "index2": index2,
            }
        )
    return samples


def read_sample_counts(input_path: pathlib.Path) -> dict[tuple[str, str], int] | None:
    """Reads per (lane, sample) from the run's only ``Demultiplex_Stats.csv``."""
    found = sorted(input_path.glob("**/Reports/Demultiplex_Stats.csv"))
    if len(found) != 1:
        logging.info(
            f"Barcode findings: {len(found)} Demultiplex_Stats.csv files found "
            "(need exactly one); counts not thresholded, putative samples skipped."
        )
        return None
    with open(found[0], newline="") as handle:
        return {
            (_lane(row["Lane"]), row["SampleID"]): int(row["# Reads"])
            for row in csv.DictReader(handle)
        }


def _finding(lane, label, index, count, sample=None, other=None, matched=""):
    return {
        "Lane": lane,
        "Label": label,
        "Undetermined_Index": index,
        "Count": count,
        "Sample_ID": sample["Sample_ID"] if sample else None,
        "Other_Sample_ID": other["Sample_ID"] if other else None,
        "Sample_Project": sample["Sample_Project"] if sample else None,
        "Matched": matched,
    }


def _explain(lane, index, count, samples) -> list[dict]:
    """Index Mixing, Index Swap, Reverse and Complement findings of one barcode."""
    parts = index.split("+")
    findings = []
    if len(parts) == 2:
        for a in samples:
            for b in samples:
                if a is not b and a["index1"] == parts[0] and b["index2"] == parts[1]:
                    findings.append(
                        _finding(
                            lane,
                            "Index Mixing",
                            index,
                            count,
                            a,
                            b,
                            f"i7 of {a['Sample_ID']}, i5 of {b['Sample_ID']}",
                        )
                    )
        for s in samples:
            if (
                s["index2"]
                and parts[0] in (s["index2"], _rc(s["index2"]))
                and parts[1] in (s["index1"], _rc(s["index1"]))
            ):
                findings.append(
                    _finding(
                        lane,
                        "Index Swap",
                        index,
                        count,
                        s,
                        matched="i7 and i5 exchanged",
                    )
                )
    for s in samples:
        for name in ("index1", "index2"):
            sample_part = s[name]
            if not sample_part:
                continue
            for label, flipped in (
                ("Reverse", sample_part[::-1]),
                ("Complement", sample_part.translate(_COMPLEMENT)),
            ):
                # A palindrome flips to itself: that is the index as given
                if flipped == sample_part:
                    continue
                for position, part in zip(("i7", "i5"), parts):
                    if part == flipped:
                        findings.append(
                            _finding(
                                lane,
                                label,
                                index,
                                count,
                                s,
                                matched=f"{position} {part} = {label.lower()} of {name} {sample_part}",
                            )
                        )
    return findings


def classify(
    samples: list[dict],
    undetermined_rows: list[dict],
    counts: dict[tuple[str, str], int] | None,
) -> list[dict]:
    """Findings for the undetermined rows (``Lane``, ``Undetermined_Index``, ``Count``)."""
    by_lane: dict[str, list[dict]] = {}
    for s in samples:
        by_lane.setdefault(s["Lane"], []).append(s)
    findings = []
    for row in sorted(undetermined_rows, key=lambda r: (r["Lane"], -r["Count"])):
        lane, index, count = row["Lane"], row["Undetermined_Index"], row["Count"]
        lane_samples = by_lane.get(lane, [])
        if not lane_samples or any(run in index for run in DISCARD_RUNS):
            continue
        # A sample's own barcode (e.g. through a collision) is not a finding
        own = {
            s["index1"] + (f"+{s['index2']}" if s["index2"] else "")
            for s in lane_samples
        }
        if index in own:
            continue
        lane_counts = None
        if counts is not None:
            lane_counts = {
                s["Sample_ID"]: counts[(lane, s["Sample_ID"])]
                for s in lane_samples
                if (lane, s["Sample_ID"]) in counts
            }
            # Only barcodes larger than the lane's smallest sample matter
            if not lane_counts or count <= min(lane_counts.values()):
                continue
        found = _explain(lane, index, count, lane_samples)
        if not found and lane_counts:
            projects = {
                s["Sample_Project"] for s in lane_samples if s["Sample_Project"]
            }
            if len(projects) == 1:
                project_counts = list(lane_counts.values())
                if len(project_counts) >= 2:
                    median = statistics.median(project_counts)
                    if abs(count - median) <= PUTATIVE_STDEVS * statistics.stdev(
                        project_counts
                    ):
                        found.append(
                            {
                                **_finding(lane, "Putative Sample", index, count),
                                "Sample_Project": next(iter(projects)),
                            }
                        )
        findings.extend(found)
    return findings


def report(
    sample_sheet: pathlib.Path,
    input_path: pathlib.Path,
    undetermined_rows: list[dict],
    output_path: pathlib.Path,
    basename: str,
) -> None:
    """Write ``<basename>_barcode_findings.csv``."""
    findings = classify(
        read_samples(sample_sheet), undetermined_rows, read_sample_counts(input_path)
    )
    output_file = output_path / f"{basename}_barcode_findings.csv"
    with open(output_file, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(findings)
    logging.info(f"{len(findings)} barcode findings saved to {output_file}.")
