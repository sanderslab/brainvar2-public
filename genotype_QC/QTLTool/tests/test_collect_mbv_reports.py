#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COLLECTOR = ROOT / "bin" / "collect_mbv_reports.py"
REPORTS = ROOT / "tests" / "fixtures" / "reports"


class CollectMbvReportsTests(unittest.TestCase):
    def run_collector(self, output: Path) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                "python3", str(COLLECTOR),
                "--reports", str(REPORTS / "RNA_A.mbv.txt"), str(REPORTS / "RNA_B.mbv.txt"),
                "--output", str(output),
            ],
            check=True,
            text=True,
            capture_output=True,
        )

    def test_collects_every_native_candidate_row(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            output = Path(tempdir) / "mbv_all.tsv.gz"
            result = self.run_collector(output)
            self.assertIn("6 candidate rows from 2 MBV reports", result.stdout)
            with gzip.open(output, "rt", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 6)
            self.assertEqual({row["rna_alignment"] for row in rows}, {"RNA_A", "RNA_B"})
            self.assertEqual({row["candidate_donor"] for row in rows}, {"DONOR_A", "DONOR_B", "DONOR_C"})
            self.assertEqual(rows[0]["perc_het_consistent"], "0.975")

    def test_gzip_output_is_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            first = Path(tempdir) / "first.tsv.gz"
            second = Path(tempdir) / "second.tsv.gz"
            self.run_collector(first)
            self.run_collector(second)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_rejects_different_candidate_sets(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            changed = root / "RNA_CHANGED.mbv.txt"
            text = (REPORTS / "RNA_B.mbv.txt").read_text(encoding="utf-8")
            changed.write_text(text.replace("DONOR_C", "DONOR_D"), encoding="utf-8")
            result = subprocess.run(
                [
                    "python3", str(COLLECTOR),
                    "--reports", str(REPORTS / "RNA_A.mbv.txt"), str(changed),
                    "--output", str(root / "mbv_all.tsv"),
                ],
                check=False,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("candidate donor set differs", result.stderr)


if __name__ == "__main__":
    unittest.main()
