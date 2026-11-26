# Bulk RNA-seq quantification workflow

This directory contains a Snakemake workflow for Salmon transcript
quantification, tximport aggregation, and HTSeq gene counting. The workflow
configuration is in `config/config.yaml`; rules are under
`workflow/rules/`, and the top-level Snakefile is `workflow/Snakefile`.

## Workflow behavior

- Salmon quantification produces per-sample `quant.sf` files using the
  configured Salmon index.
- tximport combines Salmon results into gene-level count outputs. The
  Snakefile's default target is the tximport rule output.
- HTSeq gene counting is a separate rule based on aligned BAM files; it is not
  required by the default target unless explicitly requested or used as an
  input by a selected rule.

The repository includes conda environment definitions for Salmon and HTSeq.
Confirm the exact rule inputs and outputs in the rule files and update the
configuration paths for the execution environment.

## Execution

The workflow was configured for Snakemake with S3 storage support. A typical
dry run from this directory is:

```bash
snakemake --use-conda --cores 1 -n
```

Review the storage plugin setup, AWS credentials, S3 paths, Dockerfile, and
environment definitions before running. The older examples for Docker, AWS,
and GNU Parallel in project notes are environment-specific and are not
required Snakemake commands.

## References and reproducibility

The configured reference is GENCODE v41 with the project GRCh38 reference
bundle. Record the workflow revision, reference/index versions, configuration,
and software versions for each run. FASTQs, BAMs, indexes, and generated
counts are not distributed here.
