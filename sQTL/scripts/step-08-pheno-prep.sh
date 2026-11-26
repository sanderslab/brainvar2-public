#!/bin/bash

# Usage: ./step-09-pheno-prep.sh <input_counts.nochr.gz> <chr_remove.txt> <leafcutter_scripts_dir>

set -euo pipefail

INPUT_COUNTS=$1         # ex: leafcutter_perind.counts.nochr.gz
CHR_BLACKLIST=$2        # ex: chr_remove.txt
LEAFCUTTER_SCRIPTS=$3   # ex: /tools/leafcutter/scripts

SCRIPT=${LEAFCUTTER_SCRIPTS}/prepare_phenotype_table.py

echo "[INFO] Preparing phenotype table..."
python "${SCRIPT}" "${INPUT_COUNTS}" \
    --ChromosomeBlackList "${CHR_BLACKLIST}"

echo "[INFO] Done. Outputs:"
ls -1 "${INPUT_COUNTS}".phen_chr* "${INPUT_COUNTS}".qqnorm_chr* "${INPUT_COUNTS}".*