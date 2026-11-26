#!/bin/bash

# Narjes Rohani March 17 2025
# fastqc and multiqc
# c5.4xlarge
# qc for bulk rna-seq data of BrainVar2

# Set up error handling
trap 'echo "An error occurred on line $LINENO. Exiting..."; exit 1' ERR

# Enable strict error checking
set -e

# Define the input directory where the fastq files are located
parent_dir="/mnt/s3/BrainVar/"
output_dir="./bulkRNA_QC/FASTQ_output_brainvar1"  # Change this to the desired output directory

# Check if parent directory exists
if [ ! -d "$parent_dir" ]; then
  echo "Error: Parent directory $parent_dir does not exist."
  exit 1
fi

# Step 1: Run FastQC for each sample and save reports in the output directory
echo "Starting FastQC analysis..."

# Create the output directory for FastQC reports
mkdir -p "$output_dir/fastqc_reports"

# Iterate over all fastq.gz files in the parent directory
for fastq_file in "$parent_dir"/*.R1*.fastq.gz; do
  # Extract sample name from file (assuming the sample name is the part before "_R1")
  sample_name=$(basename "$fastq_file" | cut -d'.' -f1)
  echo "Processing sample: $sample_name"

  # Create output folder for the sample inside the output directory
  sample_output_dir="$output_dir/fastqc_reports/$sample_name"
  mkdir -p "$sample_output_dir"

  # Check if corresponding R2 file exists
  r2_file="${fastq_file/.R1/.R2}"
  if [ ! -f "$r2_file" ]; then
    echo "Warning: R2 fastq file for sample $sample_name not found. Skipping..."
    continue
  fi

  # Run FastQC on the R1 and R2 files and save the reports in the sample output folder
  echo "Running FastQC for $sample_name..."
  fastqc "$fastq_file" "$r2_file" -t 8 --outdir="$sample_output_dir"
done

echo "FastQC analysis completed."

# Step 2: Run MultiQC on all FastQC reports and save the summary report
echo "Running MultiQC analysis..."

multiqc "$output_dir/fastqc_reports/" -o "$output_dir/multiqc_report"

echo "MultiQC analysis completed. Reports are saved in $output_dir/multiqc_report."
