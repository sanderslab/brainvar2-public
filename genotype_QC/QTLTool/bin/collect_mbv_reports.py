#!/usr/bin/env python3
"""Collect native QTLtools MBV reports into one lossless candidate table."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import sys
from pathlib import Path
from typing import TextIO


REPORT_SUFFIX = ".mbv.txt"
REQUIRED_COLUMNS = {
    "SampleID",
    "n_het_covered",
    "n_hom_covered",
    "perc_het_consistent",
    "perc_hom_consistent",
}


class ReportError(ValueError):
    pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reports", required=True, nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def open_output(path: Path) -> TextIO:
    raw = path.open("wb")
    if path.suffix == ".gz":
        compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
        return io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    return io.TextIOWrapper(raw, encoding="utf-8", newline="")


def read_report(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.name.endswith(REPORT_SUFFIX):
        raise ReportError(f"report does not end in {REPORT_SUFFIX}: {path.name}")
    with path.open(encoding="utf-8") as handle:
        header = handle.readline().split()
        if not header:
            raise ReportError(f"empty report: {path.name}")
        missing = REQUIRED_COLUMNS.difference(header)
        if missing:
            raise ReportError(
                f"{path.name} lacks required columns: {','.join(sorted(missing))}"
            )
        rows: list[dict[str, str]] = []
        candidate_ids: set[str] = set()
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            values = line.split()
            if len(values) != len(header):
                raise ReportError(f"{path.name}:{line_number} has the wrong field count")
            row = dict(zip(header, values))
            candidate = row["SampleID"]
            if candidate in candidate_ids:
                raise ReportError(f"{path.name} contains a duplicate candidate donor")
            candidate_ids.add(candidate)
            rows.append(row)
    if not rows:
        raise ReportError(f"report has no candidate rows: {path.name}")
    return header, rows


def collect(report_paths: list[Path], output: Path) -> tuple[int, int]:
    paths = sorted((path.resolve() for path in report_paths), key=lambda path: path.name)
    if len(paths) != len({path.name for path in paths}):
        raise ReportError("report basenames are not unique")

    output.parent.mkdir(parents=True, exist_ok=True)
    expected_header: list[str] | None = None
    expected_candidates: set[str] | None = None
    total_rows = 0
    with open_output(output) as handle:
        writer = None
        for path in paths:
            header, rows = read_report(path)
            if expected_header is None:
                expected_header = header
                output_columns = ["rna_alignment", "candidate_donor"] + [
                    column for column in header if column != "SampleID"
                ]
                writer = csv.DictWriter(
                    handle, fieldnames=output_columns, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
            elif header != expected_header:
                raise ReportError(f"native report columns differ in {path.name}")

            candidates = {row["SampleID"] for row in rows}
            if expected_candidates is None:
                expected_candidates = candidates
            elif candidates != expected_candidates:
                raise ReportError(f"candidate donor set differs in {path.name}")

            assert writer is not None
            rna_alignment = path.name[: -len(REPORT_SUFFIX)]
            for row in rows:
                output_row = {
                    "rna_alignment": rna_alignment,
                    "candidate_donor": row["SampleID"],
                    **{column: row[column] for column in header if column != "SampleID"},
                }
                writer.writerow(output_row)
                total_rows += 1
    return len(paths), total_rows


def main() -> int:
    args = parse_args()
    try:
        report_count, row_count = collect(args.reports, args.output)
    except (OSError, ReportError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2
    print(f"Collected {row_count} candidate rows from {report_count} MBV reports.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
