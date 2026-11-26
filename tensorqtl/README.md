# TensorQTL Nextflow workflow

This directory contains a DSL2 Nextflow workflow for TensorQTL cis-QTL
mapping, input preparation, and PLINK conversion. It is a separate workflow
from `eQTL/eQTL_pipeline/`; review both implementations before choosing one.

## Entry points and modes

Run `main.nf` with `--mode qtl`, `--mode plink`, or `--mode test`. The QTL
workflow reads a sample sheet specified by `--sample_file`; its fields define
the genotype, phenotype, covariate, and output paths. `workflows/` and
`modules/` define the processes and expected input structure. The checked-in
configuration files contain environment-specific settings and should be
reviewed before execution.

Example command (run from this directory):

```bash
nextflow run main.nf -c nextflow.config \
  --mode qtl \
  --sample_file /path/to/sample_sheet.csv \
  -resume
```

Use the config and sample sheet appropriate to the analysis; this example is
not a complete runnable setup on its own. QTL runs require genotype data in the
format expected by the workflow, phenotype BED files, and covariates with
matching sample identifiers. The `plink` mode invokes PLINK2; confirm the
configured VCF and sample metadata before use.

## Repository contents

- `main.nf`, `workflows/`, and `modules/` contain the workflow implementation.
- `bin/` contains helpers for filling missing covariate values, transposing
  tables, and sorting BED files.
- `conf/` and the root configuration files contain execution settings.
- Files ending in `.bak` or `.wip` are retained development copies, not the
  active workflow entry points.

