#!/usr/bin/env bash
# ------------------------------------------------------------------
# Step 5b: Assemble a merged dataset for downstream analysis:
#   - the 10 unaffected groups, copied as-is from your original output
#   - the 4 corrected groups (Micro_prenatal, SMC_postnatal,
#     Immune_postnatal, VLMC_postnatal), copied from the rescaled-k rerun
#
# Usage:
#   ./05_assemble_corrected_dataset.sh \
#       <original_dir> <corrected_dir> <merged_outdir>
#
# Example:
#   ./05_assemble_corrected_dataset.sh \
#       ~/v3_coacc_cicero ./cicero_corrected ~/v3_coacc_cicero_corrected
# ------------------------------------------------------------------
set -euo pipefail

ORIG_DIR="${1:?Usage: $0 <original_dir> <corrected_dir> <merged_outdir>}"
CORRECTED_DIR="${2:?Usage: $0 <original_dir> <corrected_dir> <merged_outdir>}"
OUT_DIR="${3:?Usage: $0 <original_dir> <corrected_dir> <merged_outdir>}"

CORRECTED_GROUPS=(Micro_prenatal SMC_postnatal Immune_postnatal VLMC_postnatal)

mkdir -p "$OUT_DIR"

echo "Copying all original connections/ccans files from $ORIG_DIR ..."
cp "$ORIG_DIR"/cicero_*_connections.tsv "$OUT_DIR"/
cp "$ORIG_DIR"/cicero_*_ccans.tsv "$OUT_DIR"/

echo "Overwriting corrected groups from $CORRECTED_DIR ..."
for g in "${CORRECTED_GROUPS[@]}"; do
  if [[ -f "$CORRECTED_DIR/cicero_${g}_connections.tsv" && -f "$CORRECTED_DIR/cicero_${g}_ccans.tsv" ]]; then
    cp "$CORRECTED_DIR/cicero_${g}_connections.tsv" "$OUT_DIR/"
    cp "$CORRECTED_DIR/cicero_${g}_ccans.tsv" "$OUT_DIR/"
    echo "  replaced: $g"
  else
    echo "  !! WARNING: corrected files for $g not found in $CORRECTED_DIR -- keeping original (uncorrected) version"
  fi
done

echo ""
echo "Done. Merged dataset assembled at: $OUT_DIR"
echo "Point --indir at this directory for step1/step2/step3 R scripts from here on."
