#!/bin/bash

# Narjes Rohani March 17 2025
# fastqc and multiqc
# c5.4xlarge
# qc for bulk rna-seq data of BrainVar2

# Set up error handling
trap 'echo "An error occurred on line $LINENO. Exiting..."; exit 1' ERR

# Enable strict error checking
set -e

# Define the input directory where the bulkRNA sample folders are located; Sestan lab directory that is mounted
parent_dir="/mnt/s3/BrainVar2/Bulk_RNA/FASTQ"
output_dir="./bulkRNA_QC/FASTQ_output"  # Change this to the desired output directory

# Check if parent directory exists
if [ ! -d "$parent_dir" ]; then
  echo "Error: Parent directory $parent_dir does not exist."
  exit 1
fi

# Step 1: Run FastQC for each sample and save reports in their respective directories
echo "Starting FastQC analysis..."

for folder in "$parent_dir"/*/; do
  # Extract sample name from folder (assuming the folder name is the sample name)
  sample_name=$(basename "$folder")
  echo "Processing sample: $sample_name"

  # Create output folder for the sample inside the output directory
  sample_output_dir="$output_dir/$sample_name"
  mkdir -p "$sample_output_dir/fastqc_reports"

  # Check if fastq files exist in the folder
  if [ ! -f "$folder"/*_R1*.fastq.gz ] || [ ! -f "$folder"/*_R2*.fastq.gz ]; then
    echo "Warning: Fastq files for sample $sample_name not found in $folder. Skipping..."
    continue
  fi

  # Run FastQC on the R1 and R2 files and save the reports in the new folder
  echo "Running FastQC for $sample_name..."
  fastqc "$folder"/*_R1*.fastq.gz "$folder"/*_R2*.fastq.gz -t 8 --outdir="$sample_output_dir/fastqc_reports"

done

echo "FastQC analysis completed."

# Step 2: Run MultiQC on all FastQC reports and save the summary report
echo "Running MultiQC analysis..."

multiqc "$output_dir"/*/fastqc_reports/ -o "$output_dir/multiqc_report"

echo "MultiQC analysis completed. Reports are saved in $output_dir/multiqc_report."
