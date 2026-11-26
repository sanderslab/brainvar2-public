#!/usr/bin/env bash
# ------------------------------------------------------------------
# Print the three files needed to evaluate whether the k-correction
# (step 5) worked, after running steps 3-6 of the corrected-dataset
# pipeline:
#   3. step1_cicero_qc.R           -> cicero_qc_plots_corrected/
#   4. approach1_cellcount_regression.R -> approach1_regression_check_corrected/
#   5. step2_cicero_crosscelltype.R
#   6. step3_cicero_prenatal_postnatal.R -> cicero_dev_comparison_corrected/
#
# Usage:
#   ./print_corrected_results.sh
#   (edit the paths below if your --outdir values differ from the
#    step-by-step plan)
# ------------------------------------------------------------------
set -euo pipefail

QC_DIR="cicero_qc_plots_corrected"
REGRESSION_DIR="approach1_regression_check_corrected"
DEV_DIR="cicero_dev_comparison_corrected"

print_file() {
  local label="$1"
  local path="$2"
  echo "================================================================"
  echo "$label"
  echo "  ($path)"
  echo "================================================================"
  if [[ -f "$path" ]]; then
    column -t -s $'\t' "$path"
  else
    echo "  !! FILE NOT FOUND -- check the step that should have produced it ran successfully"
  fi
  echo ""
}

print_file "1) Corrected step1 summary stats (all 14 groups)" \
  "$QC_DIR/step1_summary_stats.tsv"

print_file "2) Cell-count regression residuals (corrected data)" \
  "$REGRESSION_DIR/residuals_summary.tsv"

print_file "3) Micro lineage prenatal vs postnatal summary" \
  "$DEV_DIR/ALL_lineages_dev_summary.tsv"

print_file "4) Micro edge overlap (shared / prenatal-only / postnatal-only)" \
  "$DEV_DIR/Micro_edge_overlap_summary.tsv"

echo "================================================================"
echo "Done. Copy/paste all of the above."
echo "================================================================"
