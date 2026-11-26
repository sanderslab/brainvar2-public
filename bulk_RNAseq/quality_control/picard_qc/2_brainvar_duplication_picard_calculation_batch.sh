#!/bin/bash

# April 2, 2025 by Narjes Rohani

# This script calculates duplication rate of each sample of brainvar1/2 with ability to submit multiple jobs.
# Path to the directory containing BAM files with group header calculated by 1_brainvar_add_read_group_bams.sh
bam_dir="/mnt/s3/Sanders_lab/BrainVar2/bulkRNA_qc/brainvar1/output_bams_group_header/"
out_dir="./"
# Create or empty the output file
echo "File Name - Duplication Rate (%)" > $out_dir/duplication_rates.txt

# Maximum number of concurrent jobs
max_jobs=5
current_jobs=0

# Loop through all BAM files in the specified directory
for bam_file in $bam_dir/*.bam; do
  # Get the base name of the BAM file (without the path and extension)
  base_name=$(basename "$bam_file" .bam)

  # Check if the duplication metrics file already exists
  metrics_file="${out_dir}/${base_name}_duplication_metrics.txt"
  if [ -f "$metrics_file" ]; then
    echo "Metrics file for $base_name already exists. Skipping this BAM file."
    continue
  fi

  # Run Picard's MarkDuplicates for each BAM file in the background using nohup
  nohup picard MarkDuplicates \
    I="$bam_file" \
    O="${out_dir}/${base_name}_dedup.bam" \
    M="$metrics_file" \
    REMOVE_DUPLICATES=false &

  # Increment the job counter
  current_jobs=$((current_jobs + 1))
  
  # If the number of running jobs reaches the limit, wait for one to finish before continuing
  if [ "$current_jobs" -ge "$max_jobs" ]; then
    # Wait for any background job to finish
    wait -n
    # Decrement the job counter
    current_jobs=$((current_jobs - 1))
  fi

  # Store the process ID for each nohup job (optional)
  echo "Started job for $base_name with PID $!" >> $out_dir/jobs_running.txt

done

# Wait for all remaining background jobs to finish
wait

# Extract duplication rates and append them to the output file
for bam_file in $bam_dir/*.bam; do
  base_name=$(basename "$bam_file" .bam)

  # Extract the duplication rate from the metrics file and append to the summary file
  metrics_file="${out_dir}/${base_name}_duplication_metrics.txt"
  if [ -f "$metrics_file" ]; then
    duplication_rate=$(grep "Duplicate Rate" "$metrics_file" | awk '{print $3}')
    echo "$base_name - $duplication_rate" >> $out_dir/duplication_rates.txt
  else
    echo "Metrics file for $base_name not found. Skipping rate extraction." >> $out_dir/duplication_rates.txt
  fi
done

echo "Duplication rates have been saved to $out_dir/duplication_rates.txt."
