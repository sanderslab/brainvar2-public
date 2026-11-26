#!/usr/bin/env bash

# Exit on error and undefined variables
set -euo pipefail

# ==================================
# Usage example:
# nohup bash annotate_intron.sh \
#   ~/outputs/gencode_v41 \
#   ~/outputs/leafcutter_perind.counts.gz \
#   ~/outputs > annotate_intron.log 2>&1 &
#
# Arguments:
#   1: Annotation prefix (e.g., ~/outputs/gencode_v41)
#   2: Counts file (e.g., ~/outputs/leafcutter_perind.counts.gz)
#   3: Output directory (e.g., ~/outputs)
# ==================================

ANNOT_PREFIX=$1     # Prefix of annotation files
COUNTS_FILE=$2      # Leafcutter counts file
OUTDIR=$3           # Where to save leafviz.RData

# Optional: Load R module if using environment modules
# module load R/3.6.0

# Run annotate_intron.R
Rscript annotate_intron.R \
  --annot "${ANNOT_PREFIX}" \
  --counts "${COUNTS_FILE}" \
  --outdir "${OUTDIR}/" \
  --ncores 8

echo "Annotation completed. Output saved to ${OUTDIR}/leafviz.RData"