# Bulk RNA-seq alignment workflow

This directory contains a DSL2 Nextflow workflow for STAR alignment. The
entry point is `main.nf`, which runs the `align` mode and includes the
implementation in `workflows/alignment.nf`. Default execution settings,
including AWS Batch, work storage, STAR index, GENCODE annotation, and thread
count, are in `nextflow.config` and `conf/`.

## Running

Run from this directory with a sample sheet matching the fields consumed by
the workflow and with AWS Batch access configured:

```bash
nextflow run main.nf \
  --mode align \
  --sample_file /path/to/samples.csv \
  -resume
```

Review `workflows/alignment.nf`, `modules/star.nf`, and the configuration to
confirm the required sample-sheet columns, input FASTQ locations, output
layout, reference files, and resource settings for the run. The checked-in
configuration uses project S3 locations and an AWS Batch role; it is not a
local execution configuration.

## Reproducibility

Record the repository revision, Nextflow and STAR versions, reference/index
release, sample sheet, and effective configuration with each run. Controlled
FASTQs and generated alignments are not included in this repository.
