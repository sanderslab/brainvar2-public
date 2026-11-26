# Differential accessibility analysis

This folder contains the scripts used for ArchR-based differential accessibility analyses of the BrainVar single-nucleus ATAC-seq data. RNA-derived cell-type annotations are carried into the ATAC projects. Marker analyses use the ArchR `PeakMatrix`, Wilcoxon testing, bias correction for `TSSEnrichment` and `log10(nFrags)`, and a maximum of 10,000 cells per group where configured.

## DAR analysis scripts

Run the scripts with the corresponding ArchR projects and metadata. The paths in the scripts refer to the original analysis environment.

1. [`02_DAR_analyses/01_celltype_markers_by_stage.R`](02_DAR_analyses/01_celltype_markers_by_stage.R) calculates marker peaks for broad, major, and fine cell-type annotations separately in prenatal and postnatal projects. Fine cell types are selected using a minimum of 50 cells.
2. [`02_DAR_analyses/02_neuron_vs_non_neuron_markers.R`](02_DAR_analyses/02_neuron_vs_non_neuron_markers.R) groups neuronal populations and evaluates marker peaks for neuronal and non-neuronal cell groups.
3. [`02_DAR_analyses/03_prenatal_vs_postnatal_by_celltype.R`](02_DAR_analyses/03_prenatal_vs_postnatal_by_celltype.R) combines prenatal and postnatal ArchR projects, assigns developmental groups from period labels, creates group coverages, calls reproducible peaks with MACS2, builds a PeakMatrix, and evaluates prenatal/postnatal marker peaks within broad and fine cell types.
4. [`02_DAR_analyses/04_celltype_and_development_markers.R`](02_DAR_analyses/04_celltype_and_development_markers.R) creates combined cell-type and developmental-stage labels (`celltype_dev`) for ExN, InN, Micro, and OPC populations and extracts marker results using ArchR.

## Supporting scripts

- [`01_preprocessing/01_prepare_prenatal_atac.R`](01_preprocessing/01_prepare_prenatal_atac.R) and [`01_preprocessing/02_prepare_postnatal_atac.R`](01_preprocessing/02_prepare_postnatal_atac.R) prepare prenatal and postnatal ATAC inputs.

Annotation and peak-support steps described in earlier project notes are not
present in this checkout; no scripts for those steps are listed here.

## Developmental grouping

For prenatal/postnatal comparisons, periods 3–6 are designated prenatal, period 7 is transitional, and periods 8 and later are designated postnatal. The developmental comparison scripts use the available cell-type labels and stage metadata in the ArchR projects.

## Software

The primary analysis package is ArchR (v1.x), with MACS2 used for reproducible peak calling. The scripts also use R packages including `dplyr`, `data.table`, `future`, and the hg38 BSgenome reference as required by individual workflows.
