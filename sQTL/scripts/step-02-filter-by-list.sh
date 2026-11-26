#!/bin/bash

# Usage: ./step-02-filter-by-list.sh <found_sj_list> <output_file>

set -euo pipefail

SJ_LIST=$1         # e.g., found_sj_files.txt
OUTPUT_FILE=$2     # e.g., filtered_junctions.tsv

echo "[INFO] Reading SJ.out.tab files from list: $SJ_LIST"

cat "$SJ_LIST" \
| xargs cat \
| awk '($5 > 0 && $7 > 2 && $6==0 && $1!="chrM" ){out[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]++} END{for(var in out){print var}}' \
| sort -t $'\t' -k1,1 -k2,2n -k3,3n \
> "$OUTPUT_FILE"

echo "[INFO] Filtered junctions written to: ${OUTPUT_FILE}"