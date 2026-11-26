#!/bin/bash

# Usage: ./step-06-leafcutter-cluster.sh <juncfiles.txt> <outdir> <leafcutter_path>

set -euo pipefail

JUNCFILE_LIST=$1         # ex: juncfiles.txt
OUTDIR=$2                # ex: ./results/cluster
LEAFCUTTER_DIR=$3        # ex: /tools/leafcutter/

# locates leafcutter_cluster.py 
CLUSTER_SCRIPT=${LEAFCUTTER_DIR}/clustering/leafcutter_cluster.py

# Environment setting (Python 2.7)
# source activate leafcutter_env  # in case needed

mkdir -p "${OUTDIR}"

echo "[INFO] Running LeafCutter clustering..."
python "${CLUSTER_SCRIPT}" \
    -j "${JUNCFILE_LIST}" \
    -m 50 \
    -l 500000 \
    -r "${OUTDIR}"

echo "[INFO] Done. Output files:"
ls -lh "${OUTDIR}"/leafcutter_*