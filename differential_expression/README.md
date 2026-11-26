# Differential expression analysis

This folder contains the Seurat scripts for RNA differential-expression and marker analyses described in the Methods. The files are grouped by workflow stage; the numeric prefixes indicate their intended sequence.

## Workflow

### 1. Prepare data

- `01_preparation/01_downsampling_prenatal_postnatal.R` loads period-specific snRNA-seq data, applies the cell QC filters, groups cells by fine cell type, and samples up to 3,000 cells per cell type for computationally intensive downstream analyses (both downsampled and whole dataset were used; both versions when computing system is not available). Its current period setting is prenatal (periods 3–6); configure and run the postnatal counterpart as needed for periods 8 and above.

### 2. Run DEG analyses

- `02_DEG_analyses/02_celltype_markers_by_stage.R` identifies broad, major-class, and fine cell-type markers separately in prenatal and postnatal datasets.
- `02_DEG_analyses/03_neuron_vs_non_neuron_markers.R` compares neuronal and non-neuronal cell groups separately by developmental stage.
- `02_DEG_analyses/04_prenatal_vs_postnatal_by_celltype.R` compares prenatal and postnatal cells within broad cell types and fine cell subtypes, testing groups represented at both stages.
- `02_DEG_analyses/05_celltype_by_development_markers.R` identifies markers for the combined cell-type and developmental-stage groups. The intended groups are ExN, InN, Glia, and NNC; period 7 is retained as a transition group in the current script.

