# Bulk RNA-seq duplication metrics

The shell scripts in this directory add or replace BAM read groups, run
Picard MarkDuplicates to generate duplication metrics, and combine the
per-sample values into a summary. The batch and single-process scripts are
alternative execution modes; choose one according to the available compute
environment.

## Requirements and inputs

The workflow requires coordinate-sorted BAM files, samtools, Picard, and a
Unix-like shell. The scripts contain project-specific input and output paths;
review them and configure authorized paths before running. Ensure each BAM
has the intended sample/read-group metadata and that temporary and metrics
directories have sufficient space.

## Steps

1. Run `1_brainvar_add_read_group_bams.sh` if read groups must be added.
2. Run either `2_brainvar_duplication_picard_calculation.sh` or
   `2_brainvar_duplication_picard_calculation_batch.sh` on the prepared BAMs.
3. Run `3_brainvar_duplication_rate_summary.sh` on the generated Picard
   metrics files.

Check each script's arguments, concurrency limits, and output filenames
before execution. BAMs and metrics are not distributed with this repository.
