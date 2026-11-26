# find ~/brainvar/01-aligned/STAR2/ -name "Aligned.sortedByCoord.out.bam" > bam_list.250703.txt

#!/bin/bash

# Usage: nohup bash step-04-wasp.sh /home/ec2-user/data/BrainVar2/Bulk_RNA/Scripts/bam_list.250703.txt /home/ec2-user/data/AlignedBam1 > wasp_run.log 2>&1 &

set -euo pipefail

INPUT_LIST=$1
OUTDIR=$2

THREADS=8

mkdir -p "${OUTDIR}"

while IFS= read -r IN_BAM
do
    SAMPLE_ID=$(basename $(dirname "${IN_BAM}"))
    printf "[INFO] Processing: %s\n" "${SAMPLE_ID}" >&2

    OUT_BAM="${OUTDIR}/${SAMPLE_ID}.STARAligned.sortedByCoord.WASPfiltered.out.bam"

    samtools view -@ ${THREADS} -h "${IN_BAM}" |
      grep -v "vW:i:[2-7]" |
      samtools view -@ ${THREADS} -b -o "${OUT_BAM}"

    printf "[INFO] Finished: %s\n" "${OUT_BAM}" >&2

done < "${INPUT_LIST}"