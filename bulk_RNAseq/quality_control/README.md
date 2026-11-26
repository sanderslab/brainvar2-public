# Bulk RNA-seq quality control

This directory contains FastQC/MultiQC shell scripts for BrainVar1 and
BrainVar2 FASTQ files and an R script for QC of bulk RNA-seq count data.
The R script calculates normalized expression summaries and sample-level QC
metrics/plots. It is a downstream count-matrix QC analysis, not a FASTQ
alignment or read-counting workflow.

## Inputs and outputs

The scripts expect paired-end FASTQs or a gene-count matrix, together with
library/subject metadata and gene annotation as applicable. Their current
paths are project-specific. Inspect and update those paths before execution.
Outputs include FastQC reports, MultiQC summaries, normalized count tables,
sample QC summaries, and plots as defined in the scripts.

## Software

FastQC and MultiQC are required for the shell workflows. The R QC script
depends on packages declared or loaded in its source, including edgeR,
DESeq2, limma, and plotting/annotation packages. Review the script and use a
recorded R environment rather than relying on automatic package installation.

FASTQs, count matrices, metadata, and generated reports are not included in
this repository. Confirm sample identity, annotation release, and count-file
version before interpreting QC results.
