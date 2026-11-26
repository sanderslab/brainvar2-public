#!/usr/bin/env Rscript

library(rlang)
library(tidyverse)
library(progress)

# ------------------------------------------------------------------------------
# Function: subset_region
# ------------------------------------------------------------------------------

subset_region <- function(
    gwas_path,
    eqtl_path,
    output_dir
) {
  # Ensure output directory exists
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Load input files
  gwas <- read_delim(gwas_path, show_col_types = FALSE)
  eqtl <- read_delim(eqtl_path, show_col_types = FALSE)
  
  # Detect Ensembl column (handles ENSEMBLID / ENSEMBL_ID / Ensembl_ID etc.)
  ensembl_cols <- names(eqtl)
  match_col <- ensembl_cols[grepl("^ensembl", ensembl_cols, ignore.case = TRUE)]
  
  if (length(match_col) == 0) {
    stop("No Ensembl ID column found in eQTL file: ", eqtl_path)
  }
  eqtl_ensembl_col <- match_col[1]
  
  # Unique list of genes
  genes <- eqtl %>%
    pull(eqtl_ensembl_col) %>%
    unique()
  
  message("Total unique genes to process: ", length(genes))
  
  processed_genes <- character(0)
  
  pb <- progress_bar$new(
    format = "Processing genes [:bar] :percent (:current/:total) ETA: :eta",
    total = length(genes),
    width = 80,
    clear = FALSE
  )
  
  genes %>% purrr::walk(~ {
    pb$tick()
    if (.x %in% processed_genes) return()
    
    eqtl_gene <- eqtl %>% dplyr::filter(.data[[eqtl_ensembl_col]] == .x)
    
    if (nrow(eqtl_gene) > 0) {
      eqtl_gene_name <- unique(eqtl_gene$GENESYMBOL)
      eqtl_min <- min(eqtl_gene$BP)
      eqtl_max <- max(eqtl_gene$BP)
      chr <- unique(eqtl_gene$CHR)
      
      gwas_sub <- gwas %>%
        dplyr::filter(CHR == chr, BP >= eqtl_min, BP <= eqtl_max)
      
      eqtl_file <- file.path(output_dir, paste0("chr", chr, "_", .x, "_", eqtl_gene_name, "_eqtl.txt"))
      gwas_file <- file.path(output_dir, paste0("chr", chr, "_", .x, "_", eqtl_gene_name, "_gwas.txt"))
      
      write_tsv(eqtl_gene, eqtl_file)
      write_tsv(gwas_sub, gwas_file)
      
      processed_genes <<- c(processed_genes, .x)
    }
  })
  
  message("Processing complete. Total genes processed: ", length(processed_genes))
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

gwas_base_dir   <- "/mnt/s3/GWAS_data/"
eqtl_base_dir   <- "/mnt/s3/eqtl/eQTL/subsetted_by_gene"
output_base_dir <- "/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/intersections/"

# GWAS list (as in your original script)
gwas_files <- list(
  "AOW"      = "munged_noINFOFilter_GRCh37_WALK_metaGWAS_09032023.txt.gz",
  "ASD"      = "munged_70INFOFilter_GRCh37_iPSYCH-PGC_ASD_Nov2017.txt.gz",
  "ALZ"      = "munged_noINFOFilter_subsetCols_GRCh37_GCST90027158.txt.gz",
  "ADHD" = "munged_70INFOFilter_GRCh37_ADHD2022_iPSYCH_deCODE_PGC.meta.txt.gz",
  "SCZ" = "munged_70INFOFilter_GRCh37_PGC3_SCZ_wave3.primary.autosome.public.v3.txt.gz",

  "AN"       = "munged_70INFOFilter_GRCh37_pgcAN2.2019-07.txt.gz",
  "EA"       = "munged_noINFOFilter_GRCh37_EA4_additive_excl_23andMe.txt.gz",
  "EPI"      = "munged_noINFOFilter_GRCh37_ILAE_All_Epi_11.8.14.fixed_Chr_Pos.txt.gz",
  "IQ"       = "munged_noINFOFilter_GRCh37_SavageJansen_2018_intelligence_metaanalysis.txt.gz",
  "MDD"      = "munged_70INFOFilter_GRCh37_pgc-mdd2025_no23andMe_div_v3-49-46-01.txt.gz",
  "MS"       = "munged_noINFOFilter_GRCh37_discovery_metav3.0.meta.txt.gz",
  "OCD"      = "munged_70INFOFilter_GRCh37_ocd_aug2017.txt.gz",
  "PD"       = "munged_noINFOFilter_GRCh37_GCST90275127.h.txt.gz",
  "TS"       = "munged_70INFOFilter_GRCh37_TS_Oct2018.txt.gz",
  "PAU"      = "munged_noINFOFilter_GRCh37_PAU_trans_ancestry_Aug2023.txt.gz",
  "BIP"      = "munged_70INFOFilter_GRCh37_bip2024_multianc_no23andMe.txt.gz",
  "PANIC"    = "munged_70INFOFilter_GRCh37_pgc-panic2019.tsv.gz",
  "DD"       = "munged_noINFOFilter_GRCh37_supplementary_data3_GWAS_meta.txt.gz",
  "GEPI"     = "munged_noINFOFilter_GRCh37_GCST90271609.tsv.gz",
  "CEPI"     = "munged_noINFOFilter_GRCh37_GCST90271616.tsv.gz",
  "HEIGHT"   = "munged_noINFOFilter_GRCh37_GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_ALL.txt.gz",
  "SUD"      = "munged_noINFOFilter_GRCh37_Hatoum2023AddictionEuropean.txt.gz",
  "PTSD"     = "munged_noINFOFilter_GRCh37_trans_ptsd_pcs_v4_aug3_2021.tsv.gz",
  "OCD2025"  = "munged_70INFOFilter_GRCh37_daner_OCD_full_wo23andMe_190522.txt.gz",
  "FTLD_TDP" = "munged_noINFOFilter_GRCh37_GCST90558311.h.tsv.gz",
  "ALS"      = "munged_noINFOFilter_GRCh37_GCST90027163_buildGRCh37.tsv.gz"
)

# ------------------------------------------------------------------------------
# Build eqtl_datasets dynamically
# ------------------------------------------------------------------------------

build_eqtl_datasets <- function(gwas_name) {
  eqtl_dir <- file.path(eqtl_base_dir, gwas_name)
  eqtl_files <- list.files(eqtl_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
  
  lapply(eqtl_files, function(eqtl_file) {
    list(
      name = tools::file_path_sans_ext(basename(eqtl_file)),
      eqtl_file = eqtl_file
    )
  })
}

# ------------------------------------------------------------------------------
# Safe runner
# ------------------------------------------------------------------------------

run_subset_region_safe <- function(gwas_name, gwas_file, eqtl_dataset) {
  tryCatch({
    cat(sprintf("Processing: %s x %s\n", gwas_name, eqtl_dataset$name))
    
    gwas_path <- file.path(gwas_base_dir, gwas_file)
    eqtl_path <- eqtl_dataset$eqtl_file
    output_dir <- file.path(output_base_dir, paste(gwas_name, eqtl_dataset$name, sep = "_"))
    output_dir <- file.path( output_base_dir,paste(gwas_name, sub("\\.tsv$", "", eqtl_dataset$name), sep = "_"))
    if (!file.exists(gwas_path)) {
      warning(sprintf("GWAS file not found: %s", gwas_path)); return(FALSE)
    }
    if (!file.exists(eqtl_path)) {
      warning(sprintf("eQTL file not found: %s", eqtl_path)); return(FALSE)
    }
    
    subset_region(
      gwas_path = gwas_path,
      eqtl_path = eqtl_path,
      output_dir = output_dir
    )
    
    cat(sprintf("✓ Completed: %s x %s\n", gwas_name, eqtl_dataset$name))
    return(TRUE)
    
  }, error = function(e) {
    cat(sprintf("✗ Error in %s x %s: %s\n", gwas_name, eqtl_dataset$name, e$message))
    return(FALSE)
  })
}

# ------------------------------------------------------------------------------
# Main execution
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide GWAS name as argument, e.g. Rscript run_coloc.R IQ")
}

gwas_name <- args[1]
if (!gwas_name %in% names(gwas_files)) {
  stop("GWAS name not found in gwas_files: ", gwas_name)
}

gwas_file <- gwas_files[[gwas_name]]
eqtl_datasets <- build_eqtl_datasets(gwas_name)

for (eqtl_dataset in eqtl_datasets) {
  run_subset_region_safe(gwas_name, gwas_file, eqtl_dataset)
}
