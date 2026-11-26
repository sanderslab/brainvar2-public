#!/bin/bash

# Usage: bash step-01-star-1stpass.sh <samples.txt> <output_directory>

set -euo pipefail

# ---------------------------
# Input arguments
# ---------------------------
INPUT_FILE=$1   # sample list (tab delimited: sample_id, R1, R2)
OUTDIR=$2       # output directory

# ---------------------------
# Environment settings
# ---------------------------
THREADS=12
STAR_INDEX=/home/ec2-user/ref/STAR_GRCH38.p13_GENCODEv41_GTEx_no_alt-hla-decoy-contigs
GTF_FILE=/home/ec2-user/ref/gencode.v41.annotation.gtf
STAR=/home/ec2-user/miniconda/envs/leafcutter_env/bin/STAR

# ---------------------------
# Iterate over samples
# ---------------------------
while IFS=$'\t' read -r SAMPLE_ID R1 R2
do
    echo "[INFO] Processing sample: ${SAMPLE_ID}"
    
    SAMPLEPATH=${OUTDIR}/${SAMPLE_ID}
    mkdir -p ${SAMPLEPATH}
    
    ${STAR} --genomeDir ${STAR_INDEX} \
        --runThreadN ${THREADS} \
        --readFilesIn ${R1} ${R2} \
        --readFilesCommand gunzip -c \
        --outSAMattributes NH HI AS NM MD nM \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --sjdbGTFfile ${GTF_FILE} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --outFileNamePrefix ${SAMPLEPATH}/${SAMPLE_ID}.STAR \
        --limitSjdbInsertNsj 1000000
    
    echo "[INFO] Finished sample: ${SAMPLE_ID}"
done < "${INPUT_FILE}"