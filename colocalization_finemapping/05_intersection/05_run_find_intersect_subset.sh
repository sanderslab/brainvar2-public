#!/bin/bash
set -euo pipefail

# Path to R script
R_SCRIPT="./find_intersect_eqtl_subset.R"

# Logs directory
mkdir -p logs

# GWAS list (must match run_coloc.R)
GWAS_LIST=(
  AOW ASD ALZ AN EA EPI IQ MDD MS OCD PD TS PAU BIP PANIC DD GEPI CEPI HEIGHT SUD PTSD OCD2025 FTLD_TDP ALS
)
GWAS_LIST=(
  SCZ ADHD IQ 
)
for GWAS in "${GWAS_LIST[@]}"; do
  echo "Launching ${GWAS}..."
  nohup Rscript "${R_SCRIPT}" "${GWAS}" > "logs/${GWAS}.out" 2> "logs/${GWAS}.err" &
done

echo "✅ All jobs submitted with nohup (running in background)."
echo "Check progress with: tail -f logs/<GWAS>.out"
