# Genotype quality control

This directory contains Hail scripts for cohort variant QC, rare-variant
filtering, gnomAD annotation, and ancestry PCA, plus a separate QTLtools MBV
workflow under [`QTLTool/`](QTLTool/README.md).

## Hail scripts

- `hail_comVarQC_single.py` and `hail_comVarQC_single_AF0.02.py` implement
  common-variant QC/PCA branches with different allele-frequency settings.
- `hail_rareVar_HQ.py` prepares the high-quality rare-variant set and PCA
  outputs.
- `hail_get_VarQC_gnomadAnno_circos_input.py` creates a Hail MatrixTable with
  variant annotations and prepares Circos-related output.
- `hail_pca_for_anc_perdict.py` prepares ancestry PCA inputs using the
  configured reference panel.

These scripts contain Google Cloud/Hail paths and reference-resource
locations. Review the filtering logic, genome build, resource versions,
metadata fields, and output URIs in the source before running. Use an
appropriately configured Hail/Dataproc environment; the scripts are not a
portable command-line pipeline and the cohort data are not included here.
