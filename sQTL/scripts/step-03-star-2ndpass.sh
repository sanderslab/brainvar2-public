#!/bin/bash

# Usage: ./step-03-star-2ndpass.sh <samples.txt> <output_dir> <filtered_sj.tab>

set -euo pipefail

INPUT_FILE=$1
OUTDIR=$2
FILTERED_SJ=$3

# STAR setting
THREADS=16
STAR_INDEX=/home/ec2-user/ref/STAR_GRCH38.p13_GENCODEv41_GTEx_no_alt-hla-decoy-contigs
GTF_FILE=/home/ec2-user/ref/gencode.v41.annotation.gtf
STAR=/home/ec2-user/miniconda/envs/leafcutter_env/bin/STAR

# ---------------------------
# Iterate over each sample
# ---------------------------
while IFS=$'\t' read -r SAMPLE_ID BULK_ID R1 R2 VCF_ID VCF_FILE
do
    echo "[INFO] Processing ${SAMPLE_ID}"

    SAMPLEPATH=${OUTDIR}/${SAMPLE_ID}
    mkdir -p ${SAMPLEPATH}

    ${STAR} --genomeDir ${STAR_INDEX} \
        --runThreadN ${THREADS} \
        --readFilesIn ${R1} ${R2} \
        --readFilesCommand gunzip -c \
        --outSAMattributes NH HI AS NM MD nM vW \
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
        --outFileNamePrefix ${SAMPLEPATH}/${SAMPLE_ID}.STAR2 \
        --limitSjdbInsertNsj 1000000 \
        --waspOutputMode SAMtag \
        --sjdbFileChrStartEnd ${FILTERED_SJ} \
        --varVCFfile /home/ec2-user/data/BrainVar2/WGS/${VCF_FILE}

    echo "[INFO] Finished ${SAMPLE_ID}"
done < "${INPUT_FILE}"