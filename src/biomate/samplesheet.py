"""Shared Illumina sample sheet reader."""

import csv
import io
import pathlib

import polars

DATA_SECTIONS = ("[BCLConvert_Data]", "[Data]")


def read_section(
    path: pathlib.Path,
    names: tuple[str, ...] = DATA_SECTIONS,
    has_header: bool = True,
    infer_schema: bool = True,
    bare_header: tuple[str, ...] = (),
) -> polars.DataFrame:
    """Return one section of a sample sheet as a DataFrame.

    The section starts after the first line whose first field is in ``names``
    and ends at the next ``[...]`` line. Blank lines and trailing commas are
    ignored; a row of commas only becomes a row of nulls. With ``infer_schema=False`` every column is a string.

    ``bare_header`` accepts sheets without section markers: the section then
    starts at the first line whose first field is one of these column names.

    Raises:
        ValueError: If the section is missing or has no lines.
    """
    rows: list[list[str]] = []
    found = False
    with open(path, newline="") as handle:
        for fields in csv.reader(handle):
            blank_line = not fields
            while fields and not fields[-1].strip():
                fields.pop()
            if not fields:
                # A row of commas only is kept: it still takes a sample number
                if found and rows and not blank_line:
                    rows.append(fields)
                continue
            first = fields[0].strip()
            if not found:
                if first in names:
                    found = True
                elif first in bare_header:
                    found = True
                    rows.append(fields)
                continue
            if first.startswith("["):
                break
            rows.append(fields)

    if not found:
        raise ValueError(
            f"Could not find section {' or '.join(names)} in Sample Sheet {path}."
        )
    if not rows:
        raise ValueError(f"No data rows found in section {' or '.join(names)}.")

    # Same width for every row: a trailing empty field was stripped above
    width = len(rows[0]) if has_header else max(len(row) for row in rows)
    text = io.StringIO()
    csv.writer(text).writerows((row + [""] * width)[:width] for row in rows)
    return polars.read_csv(
        io.StringIO(text.getvalue()), has_header=has_header, infer_schema=infer_schema
    )
