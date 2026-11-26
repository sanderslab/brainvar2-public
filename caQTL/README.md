# caQTL input preparation

This directory contains R scripts for preparing prenatal and postnatal
pseudobulk ATAC-seq phenotypes and covariates for cis-caQTL mapping. The
scripts are numbered by workflow stage where practical, but several filenames
include historical dates and the scripts retain project-specific paths.
Review each script's inputs and output locations before running it.

## Workflow

1. **Prepare metadata and filter peaks.** The prenatal and postnatal
   `metadata_*_organize_*.R` scripts combine library, subject, and cell
   metadata. `remove_blacklist_20250626.R` removes peaks overlapping the
   configured hg38 ENCODE blacklist.
2. **Normalize pseudobulk accessibility.** The `preprocessing_*_norm_*.R`
   scripts prepare donor-level matrices and apply the normalization steps
   implemented in code (including CPM/log-CPM and quantile normalization).
   Inspect these scripts for donor exclusions and cell-type-specific filters.
3. **Prepare genomic annotations and sample covariates.** The
   `organize_genomic_region_20250626.R` and `*_peak_locs_*.R` scripts prepare
   peak coordinates/annotations. The `prepare_caQTL_input_*_20250730.R`
   scripts assemble covariates and resolve sample identifiers using the
   configured metadata.
4. **Format TensorQTL phenotypes.** The
   `*_organize_final_input_20250723.R` scripts combine normalized peak values
   with coordinates, match samples to covariates, filter peaks with more than
   20% zero values, and write phenotype BED files. Exact behavior and outputs
   are defined by the scripts.
5. **Optional follow-up utilities.** `celltype_sample_list_20250729.R`,
   `check_all_zero_peaks_20250626.R`, and `caQTL_annotation_peaks.R` support
   sample selection, diagnostics, or annotation of association results.

## Inputs and execution

The code expects project pseudobulk `SummarizedExperiment` objects,
library/subject/cell metadata, genotype and ATAC principal components, gene
annotations, and (for blacklist removal) a compatible hg38 blacklist file.
The scripts currently use fixed project-environment paths. They are not
configured as a portable command-line pipeline; update paths and verify
expected files before execution.

The final-input scripts write TensorQTL-compatible BED phenotype files and
matching covariate tables to paths configured in the R code.
