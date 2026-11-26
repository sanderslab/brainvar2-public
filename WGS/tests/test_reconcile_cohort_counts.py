#!/usr/bin/env python3

from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SPEC = spec_from_file_location(
    "reconcile_cohort_counts", ROOT / "bin" / "reconcile_cohort_counts.py"
)
assert SPEC and SPEC.loader
MODULE = module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CohortReconciliationTest(unittest.TestCase):
    def test_aggregate_reconciliation_and_draft_comparison(self) -> None:
        fixtures = ROOT / "tests" / "fixtures" / "cohort"
        report = MODULE.build_report(
            fixtures / "joint_manifest.csv",
            fixtures / "library_metadata_v1.tsv",
            fixtures / "subject_metadata.tsv",
            fixtures / "library_metadata_v2.tsv",
        )

        libraries = report["library_metadata"]
        self.assertEqual(libraries["wgs_records"], 6)
        self.assertEqual(libraries["unique_linked_donor_ids"], 4)
        self.assertEqual(libraries["unassigned_wgs_records"], 1)
        self.assertEqual(libraries["extra_wgs_records_for_resequenced_donors"], 1)
        self.assertTrue(libraries["reconciliation"]["equation_holds"])
        self.assertEqual(
            libraries["reconciliation"][
                "putative_donor_entities_if_each_unassigned_record_is_distinct"
            ],
            5,
        )
        self.assertTrue(
            report["joint_manifest"]["complete_sample_by_contig_matrix"]
        )
        self.assertTrue(
            report["cross_checks"]["manifest_and_wgs_library_sets_identical"]
        )
        comparison = report["draft_library_metadata_comparison"]
        self.assertTrue(comparison["wgs_metadata_identical"])
        self.assertEqual(
            comparison["modality_record_delta_draft_minus_released"]["bulkRNA"],
            -1,
        )


if __name__ == "__main__":
    unittest.main()
