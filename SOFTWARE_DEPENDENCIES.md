# Software Dependencies

This repository contains analysis scripts associated with the BrainVar2
manuscript.

Different analyses use different software tools, R packages, Python packages,
and computational environments.

The scripts are provided as reference code and do not constitute a single
automated pipeline.

## Bulk RNA-seq

Software used includes:

- STAR
- HTSeq
- FastQC
- MultiQC
- Picard
- samtools

R packages used include:

- edgeR
- DESeq2
- limma

## Single-cell and Multiome Analyses

Packages referenced in the scripts include:

- Seurat
- Signac
- ArchR
- monocle
- cicero

## QTL and Fine-mapping Analyses

Software and packages referenced include:

- QTLtools
- tensorQTL
- coloc
- susieR
- locuszoomr

## Functional and Regulatory Analyses

Software and packages referenced include:

- CellChat
- clusterProfiler
- fgsea
- HOMER
- ChIPseeker
- motifbreakR
- TFBSTools

## General R Packages

Examples include:

- data.table
- dplyr
- ggplot2
- tidyr
- patchwork
- ComplexHeatmap
- BiocParallel

## General Python Packages

Examples include:

- numpy
- pandas
- scipy
- scanpy
- anndata
- scikit-learn
- statsmodels
- matplotlib
- seaborn

## Additional Software

Other workflows reference tools including:

- GATK
- BWA
- bcftools
- PLINK/PLINK2
- MACS2
- Nextflow

## Software Versions and References

For the complete list of software packages, computational tools, versions,
and associated references used in this study, please refer to the
**Key Resources Table** in the accompanying BrainVar2 manuscript.

Software dependencies and versions vary across analyses. Users should
consult the individual scripts and available environment configuration
files for specific requirements.

## Software Licenses

All third-party software, R packages, Python packages, and computational
tools used or referenced in this repository remain subject to their
respective licenses, terms of use, and citation requirements.

Users must comply with the licensing conditions of each software package
and computational tool. Please consult the official software repositories
and documentation for applicable licenses and usage restrictions.

The inclusion of references to third-party software does not imply
redistribution of that software or modification of its original license.

Any third-party source code included in this repository remains subject
to its applicable license and copyright notices.

## Usage Notes

The provided scripts are intended as a reference for the computational
analyses described in the accompanying manuscript and do not constitute
an automated, ready-to-run pipeline.

Input file names, directory paths, reference datasets, software
environments, and project-specific settings must be adjusted accordingly
before running the code.

The original analyses used project-specific computational environments,
reference files, and data storage locations, which may not be included
in this repository.

Users should review the individual scripts for required inputs, outputs,
software dependencies, and configuration settings before execution.

## Reproducibility

This document provides a summary of software referenced in the repository
and does not represent a complete installation environment.

Reproducing the analyses may require access to the appropriate input
datasets, reference genomes, software versions, and computational
resources.

Users should consult the manuscript's **Key Resources Table** for
software versions and associated references.

The scripts are provided to document the computational methods used
in the manuscript and are not guaranteed to execute without modification
in a different computational environment.
