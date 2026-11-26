# Colocalization and fine-mapping workflow

This folder contains the scripts in the order used to prepare GWAS/eQTL inputs, select and intersect loci, make plots, and run colocalization. The earlier BrainVar1/Sieberts workflow is retained under the legacy subfolder in stage 07; it predates the availability of BrainVar2 eQTL data.

## Pipeline

1. **GWAS summary-statistic munging** — `01_gwas_munging/01_align_sumstats.R` uses `MungeSumstats` to format GWAS files. The script currently expects `gwas_files` to be defined before its loop and has dataset-specific paths/settings.
2. **Gene selection** — `02_gene_selection/02_find_genes_around_hits.R` selects GWAS variants with P < 5 × 10⁻⁸, creates a 2 Mb window around each, and writes nearby gene lists. The implementation requires a gene’s full annotation interval to fit inside the window, which is narrower than selecting any gene that overlaps ±1 Mb.
3. **eQTL subsetting** — `03_eqtl_subsetting/03_subset_all_eQTL_data_to_genes.sh` filters eQTL files to genes in the GWAS-derived lists by gene symbol or Ensembl ID.
4. **eQTL munging** — `04_eqtl_munging/04_munge_qtl_v2.R` standardizes FastQTL-style columns, runs `MungeSumstats`, and attempts GRCh38-to-GRCh37 conversion while retaining gene annotations. Check its input columns, sample-size handling, and output paths for the dataset in use.
5. **GWAS/eQTL intersection** — `05_intersection/05_find_intersect_eqtl_subset.R` creates per-gene eQTL and GWAS files; its shell launcher is `05_intersection/05_run_find_intersect_subset.sh`. The R script selects GWAS rows over the eQTL gene’s observed position range. Shared-variant and allele matching is handled in the colocalization stage.
6. **Locus plots** — `06_locus_plots/` contains focused GWAS/eQTL locus plotting scripts. Get and Set LDLINK_TOKEN in your environment.
7. **Colocalization and fine-mapping** — `07_colocalization/` is ordered internally by analysis: standard `coloc.abf`, LD preparation, SuSiE fine-mapping and `coloc.susie`, then focused gene analyses. Older BrainVar1/Sieberts scripts are under `90_legacy/`; the duplicate and trait-specific iterations are under `99_version_history/`.
