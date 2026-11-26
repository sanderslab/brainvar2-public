#!/usr/bin/env bash

# Exit on error and undefined variables
set -euo pipefail

# ===============================
# Usage instructions
# ===============================
# Example:
# nohup bash prepare_leafviz_annot.sh \
#   ~/ref/gencode.v41.annotation.gtf \
#   ~/outputs \
#   ~/leafcutter \
#   gencode_v41 > prepare_leafviz.log 2>&1 &

# ===============================
# Arguments
# ===============================
GTF=$1               # Full path to GTF file (e.g., ~/ref/gencode.v41.annotation.gtf)
OUTPUT_DIR=$2        # Directory to save output files (e.g., ~/outputs)
LEAFCUTTER_DIR=$3    # Path to Leafcutter repo (e.g., ~/leafcutter)
OUT_PREFIX=$4        # Prefix for output files (e.g., gencode_v41)

# ===============================
# Create output directory
# ===============================
mkdir -p "${OUTPUT_DIR}"

# ===============================
# Run gtf2leafcutter.pl
# ===============================
perl "${LEAFCUTTER_DIR}/leafviz/gtf2leafcutter.pl" \
    -o "${OUTPUT_DIR}/${OUT_PREFIX}" \
    "${GTF}"

# ===============================
# Done message
# ===============================
echo "Annotation files created in ${OUTPUT_DIR}:"
ls -lh "${OUTPUT_DIR}/${OUT_PREFIX}"_*.gz