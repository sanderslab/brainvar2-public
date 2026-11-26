#!/bin/bash

# April 2, 2025 by Narjes Rohani
# Define the input directory containing the BAM files of brainvar1/brainvar2 data: sestan-lab-share:BrainVar/eQTL_RNAseq/BAM/ , s3://sestan-lab-share/BrainVar2/Bulk_RNA/BAMs/

# The Picard tool needs the bam files to have read group header. Our bam files do not have it already. This script adds the read group to the bam files and saves the output bam files to be used by Picard in the next step.

input_dir="/mnt/s3/BrainVar/BAM"
output_dir="./output_bams"  # Specify your desired output directory

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through all BAM files in the input directory
for bam_file in "$input_dir"/*.bam; do
  # Extract the sample name from the file (you can adjust this if you want)
  sample_name=$(basename "$bam_file" .bam)
  
  # Extract the numeric portion of the sample name (assuming the numeric portion is in the beginning)
  num_part=$(echo "$sample_name" | sed 's/[^0-9]*\([0-9]*\).*/\1/')

 
  # Define the output BAM file path
  output_bam="$output_dir/${sample_name}_with_readgroup.bam"

  # Run samtools addreplacerg to add read group
  samtools addreplacerg -r "@RG\tID:RG1\tSM:${sample_name}\tPL:Illumina\tLB:Library.fa" -o "$output_bam" "$bam_file"

    echo "Read group added for $sample_name. Output saved to $output_bam"
 
  fi
done

echo "All BAM files have been processed."
