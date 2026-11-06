"""Reading the CSV tables swingersim writes: one rule for every reader (GUI, reports, tools).

Leading lines starting with '#' (after spaces) are metadata of the writers (provenance, units,
epoch) and are skipped, as are blank lines; the first remaining line is the header."""

from __future__ import annotations

import csv
from collections.abc import Iterator
from pathlib import Path


def data_lines(f) -> Iterator[str]:
    """The lines of an open text file without '#' comment lines and blank lines."""
    for ln in f:
        if ln.strip() and not ln.lstrip().startswith("#"):
            yield ln


def read_rows(path: str | Path) -> list[dict[str, str]]:
    """Rows of a CSV table as dicts of stripped strings (keys stripped; empty cells '')."""
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as f:
        return [
            {k.strip(): (v or "").strip() for k, v in r.items() if k}
            for r in csv.DictReader(data_lines(f))
        ]


def header(path: str | Path) -> list[str]:
    """Column names of a CSV table."""
    with Path(path).open(newline="", encoding="utf-8", errors="replace") as f:
        for ln in data_lines(f):
            return [c.strip() for c in next(csv.reader([ln]))]
    return []


__all__ = ["data_lines", "header", "read_rows"]
