#!/bin/bash

# April 2, 2025 by Narjes Rohani
 
# The same script as 2_brainvar_duplication_picard_calculation_batch.sh without the option for multi-job submission.

# Path to the directory containing BAM files
bam_dir="/mnt/s3/Sanders_lab/BrainVar2/bulkRNA_qc/brainvar1/output_bams_group_header/"
out_dir="./"
# Create or empty the output file
echo "File Name - Duplication Rate (%)" > $output_file

# Loop through all BAM files in the specified directory
for bam_file in $bam_dir/*.bam; do
  # Get the base name of the BAM file (without the path and extension)
  base_name=$(basename "$bam_file" .bam)

  # Run Picard's MarkDuplicates to generate the duplication metrics file
  picard MarkDuplicates \
    I="$bam_file" \
    O="${out_dir}/${base_name}_dedup.bam" \
    M="${out_dir}/${base_name}_duplication_metrics.txt" \
    REMOVE_DUPLICATES=false

done

echo "Duplication rates have been saved."
