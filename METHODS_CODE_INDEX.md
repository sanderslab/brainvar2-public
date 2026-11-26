# Methods-to-code index

This index maps the current repository folders to the analysis topics in the
manuscript Methods. Check each folder README.

| Methods topic | Repository folder | Scope in this checkout |
|---|---|---|
| Bulk RNA-seq processing | `bulk_RNAseq/` | Alignment, read counting, and QC scripts and workflow definitions. |
| Whole-genome sequencing | `WGS/` | Documented germline small-variant workflow; see its README for provenance and scope. |
| Genotype QC and sample matching | `genotype_QC/` | Hail QC/PCA scripts and the QTLtools MBV matching workflow. |
| GWAS gene analysis | `GWAS/` | FUMA/MAGMA workflow documentation and QQ plotting of MAGMA results. |
| Differential expression | `differential_expression/` | RNA celltype and development markers/DEGs. |
| Differential accessibility | `differential_accessibility/` | ATAC celltype and development markers/DARs. |
| Trait and pathway enrichment | `trait_enrichment/` | FGSEA, LDSC heritability analyses. |
| Shared biological programs | `biological_domains/` | NMF domains analysis and GO Biological Process annotation. |
| Motif enrichment | `motif_enrichment/` | HOMER peak input preparation, run commands, and summaries. |
| eQTL | `eQTL/`, `tensorqtl/` | Pseudobulk preparation and TensorQTL/Nextflow materials. |
| caQTL | `caQTL/` | Peak preprocessing, annotation, and TensorQTL input preparation. |
| sQTL | `sQTL/` | STAR-WASP/LeafCutter scripts, phenotype preparation, and plotting utilities. |
| Colocalization and fine-mapping | `colocalization_finemapping/` | GWAS/eQTL preparation, intersection, locus plotting, coloc, and SuSiE scripts. |
| GO functional categories | `gene_function_category/` | GO annotation/clustering, functional classification, and expression plots. |
| Cell-cell communication | `cell_cell_communication/` | Donor-level CellChat analyses and supporting diagnostics. |
| Gene regulatory networks | `gene_regulatory_network/` | Phase and locus-centered GRN inference, visualization, and variant motif analysis. |
| Chromatin co-accessibility | `CICERO/` | CICERO workflow and supporting analyses. |
