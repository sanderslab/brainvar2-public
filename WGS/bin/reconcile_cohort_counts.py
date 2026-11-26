#!/usr/bin/env python3
"""Reconcile aggregate WGS library, donor, and joint-manifest counts."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any


def read_table(path: Path) -> list[dict[str, str]]:
    delimiter = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"{path}: missing header")
        return [
            {key: (value or "").strip() for key, value in row.items()}
            for row in reader
        ]


def require_columns(rows: list[dict[str, str]], columns: set[str], path: Path) -> None:
    if not rows:
        raise ValueError(f"{path}: no data rows")
    missing = columns - rows[0].keys()
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")


def source_metadata(path: Path, row_count: int) -> dict[str, Any]:
    stat = path.stat()
    modified = datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc)
    return {
        "path": str(path.resolve()),
        "bytes": stat.st_size,
        "modified_utc": modified.isoformat(),
        "data_rows": row_count,
    }


def flag_counts(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    return dict(sorted(Counter(row[column] or "<blank>" for row in rows).items()))


def library_summary(rows: list[dict[str, str]], path: Path) -> tuple[dict[str, Any], set[str]]:
    required = {
        "LibraryID",
        "LibraryModality",
        "SubjectID",
        "ConfirmedID",
        "Usable",
    }
    require_columns(rows, required, path)
    wgs = [row for row in rows if row["LibraryModality"].upper() == "WGS"]
    library_ids = [row["LibraryID"] for row in wgs]
    if any(not item for item in library_ids):
        raise ValueError(f"{path}: blank LibraryID in a WGS row")
    if len(library_ids) != len(set(library_ids)):
        raise ValueError(f"{path}: duplicate WGS LibraryID values")

    linked_subjects = [row["SubjectID"] for row in wgs if row["SubjectID"]]
    multiplicity = Counter(linked_subjects)
    extra_records = sum(count - 1 for count in multiplicity.values())
    unassigned = len(wgs) - len(linked_subjects)
    unique_linked = len(multiplicity)

    summary = {
        "all_library_records": len(rows),
        "modality_records": dict(
            sorted(Counter(row["LibraryModality"] or "<blank>" for row in rows).items())
        ),
        "wgs_records": len(wgs),
        "unique_wgs_library_ids": len(library_ids),
        "linked_wgs_records": len(linked_subjects),
        "unique_linked_donor_ids": unique_linked,
        "unassigned_wgs_records": unassigned,
        "donors_with_multiple_wgs_libraries": sum(
            count > 1 for count in multiplicity.values()
        ),
        "extra_wgs_records_for_resequenced_donors": extra_records,
        "donors_by_wgs_library_count": dict(
            sorted(Counter(multiplicity.values()).items())
        ),
        "confirmed_id_counts": flag_counts(wgs, "ConfirmedID"),
        "usable_counts": flag_counts(wgs, "Usable"),
        "reconciliation": {
            "wgs_records": len(wgs),
            "unique_linked_donors": unique_linked,
            "unassigned_wgs_records": unassigned,
            "extra_resequencing_records": extra_records,
            "equation_holds": len(wgs) == unique_linked + unassigned + extra_records,
            "putative_donor_entities_if_each_unassigned_record_is_distinct": (
                unique_linked + unassigned
            ),
        },
    }
    return summary, set(library_ids)


def subject_summary(rows: list[dict[str, str]], path: Path) -> dict[str, Any]:
    required = {
        "SubjectID",
        "matching_WGS_library",
        "ConfirmedID_WGS",
        "Usable_DNA",
    }
    require_columns(rows, required, path)
    subject_ids = [row["SubjectID"] for row in rows]
    if any(not item for item in subject_ids):
        raise ValueError(f"{path}: blank SubjectID")
    if len(subject_ids) != len(set(subject_ids)):
        raise ValueError(f"{path}: duplicate SubjectID values")
    return {
        "subject_records": len(rows),
        "unique_subject_ids": len(subject_ids),
        "subjects_with_matching_wgs_library": sum(
            bool(row["matching_WGS_library"]) for row in rows
        ),
        "confirmed_id_wgs_counts": flag_counts(rows, "ConfirmedID_WGS"),
        "usable_dna_counts": flag_counts(rows, "Usable_DNA"),
    }


def manifest_summary(
    rows: list[dict[str, str]], path: Path
) -> tuple[dict[str, Any], set[str]]:
    require_columns(rows, {"sample", "contig"}, path)
    samples = [row["sample"] for row in rows]
    contigs = sorted({row["contig"] for row in rows})
    if any(not item for item in samples) or any(not item for item in contigs):
        raise ValueError(f"{path}: blank sample or contig")
    pair_counts = Counter((row["sample"], row["contig"]) for row in rows)
    unique_samples = set(samples)
    expected = len(unique_samples) * len(contigs)
    complete = len(rows) == expected and all(count == 1 for count in pair_counts.values())
    rows_per_sample = Counter(samples)
    return {
        "records": len(rows),
        "unique_sample_ids": len(unique_samples),
        "unique_contigs": len(contigs),
        "expected_complete_matrix_records": expected,
        "complete_sample_by_contig_matrix": complete,
        "minimum_rows_per_sample": min(rows_per_sample.values()),
        "maximum_rows_per_sample": max(rows_per_sample.values()),
    }, unique_samples


def draft_comparison(
    released_rows: list[dict[str, str]],
    draft_rows: list[dict[str, str]],
    draft_path: Path,
) -> dict[str, Any]:
    required = {"LibraryID", "LibraryModality"}
    require_columns(draft_rows, required, draft_path)
    released = {row["LibraryID"]: row for row in released_rows}
    draft = {row["LibraryID"]: row for row in draft_rows}
    if len(released) != len(released_rows) or len(draft) != len(draft_rows):
        raise ValueError("library metadata comparison requires unique LibraryID values")

    released_wgs = {
        key: row for key, row in released.items() if row["LibraryModality"].upper() == "WGS"
    }
    draft_wgs = {
        key: row for key, row in draft.items() if row["LibraryModality"].upper() == "WGS"
    }
    changed_columns: Counter[str] = Counter()
    changed_wgs = 0
    for key in released_wgs.keys() & draft_wgs.keys():
        changed = {
            column
            for column in released_wgs[key].keys() | draft_wgs[key].keys()
            if released_wgs[key].get(column, "") != draft_wgs[key].get(column, "")
        }
        if changed:
            changed_wgs += 1
            changed_columns.update(changed)

    released_modalities = Counter(row["LibraryModality"] for row in released_rows)
    draft_modalities = Counter(row["LibraryModality"] for row in draft_rows)
    return {
        "released_records": len(released_rows),
        "draft_records": len(draft_rows),
        "modality_record_delta_draft_minus_released": {
            modality: draft_modalities[modality] - released_modalities[modality]
            for modality in sorted(released_modalities.keys() | draft_modalities.keys())
        },
        "wgs_records_added": len(draft_wgs.keys() - released_wgs.keys()),
        "wgs_records_removed": len(released_wgs.keys() - draft_wgs.keys()),
        "wgs_records_changed": changed_wgs,
        "changed_wgs_columns": dict(sorted(changed_columns.items())),
        "wgs_metadata_identical": (
            released_wgs.keys() == draft_wgs.keys() and changed_wgs == 0
        ),
    }


def build_report(
    manifest_path: Path,
    library_path: Path,
    subject_path: Path,
    draft_library_path: Path | None = None,
) -> dict[str, Any]:
    manifest_rows = read_table(manifest_path)
    library_rows = read_table(library_path)
    subject_rows = read_table(subject_path)

    manifest, manifest_samples = manifest_summary(manifest_rows, manifest_path)
    libraries, wgs_libraries = library_summary(library_rows, library_path)
    subjects = subject_summary(subject_rows, subject_path)

    report: dict[str, Any] = {
        "sources": {
            "joint_manifest": source_metadata(manifest_path, len(manifest_rows)),
            "library_metadata": source_metadata(library_path, len(library_rows)),
            "subject_metadata": source_metadata(subject_path, len(subject_rows)),
        },
        "joint_manifest": manifest,
        "library_metadata": libraries,
        "subject_metadata": subjects,
        "cross_checks": {
            "manifest_samples_in_wgs_metadata": len(manifest_samples & wgs_libraries),
            "manifest_samples_absent_from_wgs_metadata": len(manifest_samples - wgs_libraries),
            "wgs_libraries_absent_from_manifest": len(wgs_libraries - manifest_samples),
            "manifest_and_wgs_library_sets_identical": manifest_samples == wgs_libraries,
        },
        "interpretation": {
            "sequenced_wgs_libraries": libraries["wgs_records"],
            "directly_linked_unique_donors": libraries["unique_linked_donor_ids"],
            "putative_donor_entities": libraries["reconciliation"][
                "putative_donor_entities_if_each_unassigned_record_is_distinct"
            ],
            "usable_wgs_library_records": libraries["usable_counts"].get("Yes", 0),
            "usable_unique_donors": subjects["usable_dna_counts"].get("Yes", 0),
            "putative_donor_entities_caveat": (
                "This total assumes each unassigned WGS record represents a distinct donor; "
                "the metadata do not provide SubjectID keys that prove that assumption."
            ),
        },
    }

    if draft_library_path:
        draft_rows = read_table(draft_library_path)
        report["sources"]["draft_library_metadata"] = source_metadata(
            draft_library_path, len(draft_rows)
        )
        report["draft_library_metadata_comparison"] = draft_comparison(
            library_rows, draft_rows, draft_library_path
        )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Report aggregate BrainVar2 WGS cohort counts without emitting identifiers."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--library-metadata", required=True, type=Path)
    parser.add_argument("--subject-metadata", required=True, type=Path)
    parser.add_argument("--draft-library-metadata", type=Path)
    parser.add_argument("--output", type=Path, help="JSON output; defaults to stdout")
    args = parser.parse_args()

    try:
        report = build_report(
            args.manifest,
            args.library_metadata,
            args.subject_metadata,
            args.draft_library_metadata,
        )
        rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.output:
            args.output.write_text(rendered, encoding="utf-8")
        else:
            sys.stdout.write(rendered)
    except (OSError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
