#!/usr/bin/env python3
"""Validate BrainVar2 workflow samplesheets."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path
import re
import sys


SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def nonempty(row: dict[str, str], columns: tuple[str, ...], line: int) -> None:
    for column in columns:
        if not row.get(column, "").strip():
            raise ValueError(f"line {line}: missing value in '{column}'")


def resolve_existing(value: str, sheet: Path, line: int, column: str) -> None:
    path = Path(value)
    candidates = [path] if path.is_absolute() else [Path.cwd() / path]
    if not any(candidate.exists() for candidate in candidates):
        raise ValueError(f"line {line}: '{column}' path does not exist: {value}")


def validate_parabricks(rows: list[dict[str, str]], sheet: Path) -> None:
    seen_pairs: set[tuple[str, str, str]] = set()
    for line, row in enumerate(rows, start=2):
        nonempty(row, ("sample", "fastq_1", "fastq_2"), line)
        sample = row["sample"].strip()
        if not SAMPLE_RE.fullmatch(sample):
            raise ValueError(f"line {line}: invalid sample identifier '{sample}'")
        pair = (sample, row["fastq_1"].strip(), row["fastq_2"].strip())
        if pair in seen_pairs:
            raise ValueError(f"line {line}: duplicate FASTQ pair for sample '{sample}'")
        seen_pairs.add(pair)
        resolve_existing(row["fastq_1"], sheet, line, "fastq_1")
        resolve_existing(row["fastq_2"], sheet, line, "fastq_2")


def validate_joint(rows: list[dict[str, str]], sheet: Path, contigs: list[str]) -> None:
    records: list[tuple[str, str]] = []
    for line, row in enumerate(rows, start=2):
        nonempty(row, ("sample", "contig", "gvcf", "gvcf_index"), line)
        sample = row["sample"].strip()
        contig = row["contig"].strip()
        if not SAMPLE_RE.fullmatch(sample):
            raise ValueError(f"line {line}: invalid sample identifier '{sample}'")
        if contig not in contigs:
            raise ValueError(f"line {line}: unexpected contig '{contig}'")
        resolve_existing(row["gvcf"], sheet, line, "gvcf")
        resolve_existing(row["gvcf_index"], sheet, line, "gvcf_index")
        records.append((sample, contig))

    counts = Counter(records)
    duplicates = [f"{sample}:{contig}" for (sample, contig), count in counts.items() if count != 1]
    if duplicates:
        raise ValueError("duplicate sample-contig records: " + ", ".join(sorted(duplicates)[:10]))

    samples = sorted({sample for sample, _ in records})
    missing = [f"{sample}:{contig}" for sample in samples for contig in contigs if (sample, contig) not in counts]
    if missing:
        raise ValueError("incomplete sample-contig matrix; missing: " + ", ".join(missing[:10]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("parabricks", "joint"))
    parser.add_argument("samplesheet", type=Path)
    parser.add_argument("--contigs", default="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM")
    args = parser.parse_args()

    contigs = [item.strip() for item in args.contigs.split(",") if item.strip()]
    if not contigs or len(contigs) != len(set(contigs)):
        parser.error("--contigs must be a non-empty unique list")

    try:
        with args.samplesheet.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
        if not rows:
            raise ValueError("samplesheet has no data rows")
        if args.stage == "parabricks":
            validate_parabricks(rows, args.samplesheet)
        else:
            validate_joint(rows, args.samplesheet, contigs)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print(f"OK: {args.stage} samplesheet has {len(rows)} valid rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
