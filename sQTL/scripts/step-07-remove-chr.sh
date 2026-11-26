#!/bin/bash

# Usage: ./step-07-remove-chr.sh <input_counts.gz> <output_counts.gz>

set -euo pipefail

INPUT_FILE=$1
OUTPUT_FILE=$2

echo "[INFO] Removing 'chr' prefix from: ${INPUT_FILE}"
gzcat "${INPUT_FILE}" | sed '2,$s/chr//' | gzip -c > "${OUTPUT_FILE}"
echo "[INFO] Output written to: ${OUTPUT_FILE}"