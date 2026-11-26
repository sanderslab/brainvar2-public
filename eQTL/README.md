# eQTL analysis code

This folder contains eQTL phenotype preparation, genotype preparation, and QTL mapping scripts. It has one Nextflow/TensorQTL pipeline (`eQTL_pipeline/`), a separate legacy FastQTL shell workflow, pseudobulk preparation scripts, and bMIND/cell-type-deconvolution documentation. It does not contain a separate `tensorqtl/` directory.

## Workflow overview

1. **Prepare expression phenotypes and covariates.** The pseudobulk scripts aggregate single-cell RNA counts by sample and broad cell type, apply cell-level QC, and prepare normalized expression BED files. The preparation notes document library metadata, genotype matching, and covariate construction.
2. **Prepare genotypes and covariates for the Nextflow pipeline.** The `plink` mode converts VCF inputs to PLINK files using PLINK2. The `prepare` workflow transposes covariate files. The QTL workflow sorts phenotype BED sample columns and fills missing covariate values with row medians using the included Python helpers.
3. **Run QTL mapping.** The TensorQTL workflow runs cis permutation mapping with 1,000 permutations and cis nominal mapping. The nominal process writes chromosome-specific parquet output files. 
4. **Cell-type deconvolution.** Notes for bMIND and related estimation examples are in `eQTL_preparation/celltype_deconvolution/`.

## Main code locations

- `eQTL_pipeline/main.nf` — Nextflow entry point. It selects `test`, `prepare`, `qtl`, or `plink` mode.
- `eQTL_pipeline/workflows/tensorqtl.nf` — prepares phenotype/covariate files and runs the two TensorQTL cis modes.
- `eQTL_pipeline/workflows/plink.nf` and `modules/plink.nf` — convert genotype VCF inputs to PLINK format.
- `eQTL_pipeline/modules/tensorqtl.nf` — TensorQTL processes and Python-based file-preparation processes.
- `eQTL_pipeline/eqtl_nextflow.config` — AWS Batch execution settings and default genotype/PLINK locations for the eQTL workflow.
- `eQTL_pipeline/caqtl_nextflow.config` — corresponding caQTL configuration.
- `eQTL_preparation/pseudobulk/` — count aggregation, TMM/CPM normalization, and pseudobulk phenotype preparation.
- `eQTL_preparation/celltype_deconvolution/` — bMIND and deconvolution notes.

## Running the Nextflow pipeline

Run from `eQTL_pipeline/` with Nextflow and the AWS Batch configuration available. The workflow requires a mode-specific CSV sample sheet and configured input files. For example:

```bash
nextflow run main.nf -c eqtl_nextflow.config \
  --mode qtl \
  --sample_file eqtl_covariates.csv \
  -resume
```

Other modes are `prepare`, `plink`, and `test`. Inspect `eqtl_nextflow.config`, the sample sheet, and the paths under `params` before submitting jobs; the configuration contains project-specific S3 locations and AWS Batch roles. The pipeline is configured for AWS Batch rather than a local executor by default.

The `plink` workflow uses `plink2 --vcf ... --psam ... --output-chr chrM --make-bed --split-par hg38`. The TensorQTL processes use cis mode with covariates; the permutation process requests 1,000 permutations, while the nominal process writes parquet results. Actual sample groups, PC sets, phenotypes, and covariates are determined by the CSV and files supplied to the run.


## Further notes

- [Preparation and read-in notes](eQTL_preparation/00_readin_prepare.md)
- [Bulk RNA QC notes](eQTL_preparation/01_bulk_data_qc.md)
- [bMIND/cell-type deconvolution notes](eQTL_preparation/celltype_deconvolution/README.md)
