# sQTL analysis scripts

This directory contains the BrainVar sQTL preparation scripts, supporting
LeafCutter utilities and codes. 

## Main processing sequence

1. STAR first-pass alignment: `step-01-star-1stpass.sh`.
2. Junction filtering: `step-02-filter-by-list.sh`.
3. STAR second-pass alignment and WASP tagging: `step-03-star-2ndpass.sh`.
4. WASP filtering: `step-04-wasp_v1.1.sh`.
5. BAM-to-junction conversion: `step-05-bam2junc_v1.1.sh`.
6. LeafCutter clustering: `step-06-leafcutter-cluster.sh`.
7. Chromosome-name adjustment: `step-07-remove-chr.sh`.
8. Phenotype preparation: `step-08-pheno-prep.sh`.
9. Chromosome-specific input preparation: `step-09-prepare_chr1_22.sh`.
10. Normalized phenotype concatenation: `step-10-concat-qqnorm.sh`.
11. Stage-specific phenotype and covariate preparation:
    `phenotype_covariate_process_v1.1.R`.

The scripts depend on STAR, WASP, samtools/bcftools, LeafCutter, R/Python
packages, and external reference and genotype files. 
## Supporting scripts

- `sQTL_mung.R`, `liftover_sqtl.py`, and `sqtl_concat_dask.py` support
  association-result processing.
- `intron_annotation/` contains LeafCutter intron annotation and summary
  scripts; see its README for details.
- `leafcutter/` contains the LeafCutter source and documentation used by the
  workflow.

