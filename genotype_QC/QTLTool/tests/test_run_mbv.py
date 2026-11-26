#!/usr/bin/env python3

from __future__ import annotations

import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "bin" / "run_mbv.sh"


class RunMbvTests(unittest.TestCase):
    def test_dry_run_preserves_historical_report_naming(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            root = Path(tempdir)
            bam = root / "RNA_A.Aligned.sortedByCoord.out.bam"
            vcf = root / "cohort.vcf.gz"
            for path in (bam, Path(str(bam) + ".bai"), vcf, Path(str(vcf) + ".tbi")):
                path.write_text("synthetic\n", encoding="utf-8")
            result = subprocess.run(
                [
                    "bash", str(RUNNER), "--vcf", str(vcf), "--outdir", str(root / "out"),
                    "--dry-run", str(bam),
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            self.assertIn("QTLtools mbv", result.stdout)
            self.assertIn("RNA_A.Aligned.sortedByCoord.out.bam.mbv.txt", result.stdout)
            self.assertIn("--filter-mapping-quality 15", result.stdout)


if __name__ == "__main__":
    unittest.main()
