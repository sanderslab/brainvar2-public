# eQTL TensorQTL Nextflow workflow

This directory contains a DSL2 Nextflow workflow with four modes: `test`,
`prepare`, `qtl`, and `plink`. The entry point is `main.nf`; mode-specific
processes are defined under `workflows/` and `modules/`. AWS Batch settings
and project-specific paths are supplied by the configuration files.

## Modes

- `test` runs the workflow's test process.
- `prepare` prepares covariate inputs using the included helper utilities.
- `qtl` runs cis permutation and nominal TensorQTL processes for the supplied
  phenotype, genotype, and covariate inputs.
- `plink` converts a VCF and sample metadata to PLINK files using PLINK2.

See `main.nf`, `workflows/`, `modules/`, and the selected config for required
parameters and channel contracts. The workflow is configured for AWS Batch;
execution requires an appropriate Nextflow installation, AWS access, and the
configured Batch resources.

## Example

Run from this directory after reviewing `eqtl_nextflow.config` and the sample
sheet. The exact parameter names and required files must match the checked-in
configuration and process definitions.

```bash
nextflow run main.nf -c eqtl_nextflow.config \
  --mode qtl \
  --sample_file eqtl_covariates.csv \
  -resume
```

For genotype conversion, select `--mode plink`; `prepare` and `test` are also
available. The sample command is illustrative and does not run a complete
analysis without the project inputs and cloud configuration.

## Inputs, outputs, and reproducibility

Inputs and output locations are defined by the sample sheet, workflow
parameters, and configuration. Review the Nextflow trace and process logs for
each run. Record the workflow revision, Nextflow and container/software
versions, genome build, phenotype and covariate files, sample matching, and
TensorQTL options alongside results. This repository does not include the
controlled genotype or phenotype datasets.
