#!/bin/bash

# Usage: ./step-11-concat-qqnorm.sh <qqnorm_prefix> <output_dir>

set -euo pipefail

PREFIX=$1        # 예: /home/ec2-user/data/cluster/leafcutter_perind.counts.nochr.gz.qqnorm
OUTDIR=$2

OUTFILE="${OUTDIR}/$(basename ${PREFIX})_all.gz"

echo "[INFO] Concatenating chr1-22 qqnorm files into: ${OUTFILE}"
rm -f "${OUTFILE}"  

# overwrite (>)
CHR=1
FILE="${PREFIX}_${CHR}.gz"
echo "  [chr${CHR}] -> ${FILE}"
zcat "${FILE}" | bgzip -c > "${OUTFILE}"

# 2~22 파일은 append (>>)
for CHR in {2..22}; do
    FILE="${PREFIX}_${CHR}.gz"
    echo "  [chr${CHR}] -> ${FILE}"
    zcat "${FILE}" | bgzip -c >> "${OUTFILE}"
done

tabix -p bed "${OUTFILE}"
echo "[INFO] Output indexed: ${OUTFILE}.tbi"