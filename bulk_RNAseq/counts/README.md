# Bulk RNA-seq gene counts

This directory contains scripts for per-library HTSeq counting and for
combining gene-count matrices. The scripts are `htseq.sh`,
`make_count_matrix.R`, `merge_bv_matrices.R`, and
`build_gene_matrix.R`/`build_gene_matrix_fastq.R`. Review each script for its
expected input paths, annotation, strandedness, and output format before use.

## Historical count-generation settings

The recorded HTSeq workflow used GENCODE v41 annotations and stranded BAMs.
The notes associated with that run specify BAM input, the `union` overlap
mode, and stranded counting. Confirm the actual flags in `htseq.sh` before
reproducing results; the checked-in scripts may represent different run
iterations.

Count matrices and BAMs were stored in project-controlled locations and are
not included in this repository. Do not rely on old storage paths in archived
notes. Set input and output paths to authorized locations and record the
annotation and software versions for each run.
