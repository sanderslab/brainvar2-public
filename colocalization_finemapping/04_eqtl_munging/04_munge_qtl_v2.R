# =============================================
# Load libraries
# =============================================
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("MungeSumstats")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

library(BSgenome.Hsapiens.NCBI.GRCh38)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(tidyverse)
library(data.table)
library(MungeSumstats)

# =============================================
# Preprocessing function for one file
# =============================================
preprocess_fastqtl <- function(file_path) {
  cat("Processing:", basename(file_path), "\n")
  dt <- fread(file_path)

  # Standardize column names
  setnames(dt, old = c("variant_id","pval_nominal", "slope", "slope_se", "af"),
           new = c("SNP","P", "BETA", "SE", "FRQ"), skip_absent = TRUE)

  # Estimate sample size if missing
  if ("ma_samples" %in% names(dt) && !"N" %in% names(dt)) {
    dt[, N := pmax(ma_samples * 2, 100)]
  }

  # Keep useful columns
  cols_for_munge <- intersect(names(dt),
                              c("SNP", "CHR", "BP", "A1", "A2",
                                "BETA", "SE", "P", "FRQ", "N"))
  cols_for_reference <- c("phenotype_id", "variant_id", "GeneSymbol", "EnsemblID", "GeneId")
  all_cols <- c(cols_for_munge, intersect(names(dt), cols_for_reference))
  return(dt[, ..all_cols])
}

# =============================================
# Function to merge munged data with gene info
# =============================================
merge_munged_with_gene_info <- function(munged_file, gene_info_file, output_file) {
  # Read munged data
  munged_dt <- fread(munged_file)
  
  # Read gene info
  gene_info_dt <- fread(gene_info_file)
  
  # Merge on SNP (which should be the linking variable)
  merged_dt <- merge(munged_dt, gene_info_dt, by = "SNP", all.x = TRUE)
  
  # Write merged file
  fwrite(merged_dt, output_file, sep = "\t")
  
  cat("  -> Merged file with gene info:", basename(output_file), "\n")
  return(merged_dt)
}

# =============================================
# Function to process ALL files in a path
# =============================================
munge_all_fastqtl <- function(eqtl_path, pattern = "\\.tsv$", keep_gene_info_merged = TRUE) {
  # Set up output dirs
  output_path <- file.path('/mnt/s3/eqtl/eQTL/', "munged")
  dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

  # Create logs dir inside munged
  logs_path <- file.path(output_path, "logs")
  dir.create(logs_path, showWarnings = FALSE, recursive = TRUE)

  # Find all .tsv files
  eqtl_files <- list.files(eqtl_path, pattern = pattern, full.names = TRUE)
  cat("Found", length(eqtl_files), "files to process\n\n")

  for (f in eqtl_files) {
    tryCatch({
      # --- Preprocess
      processed_dt <- preprocess_fastqtl(f)
      if (nrow(processed_dt) == 0) {
        cat("  No rows, skipping...\n")
        next
      }

      # --- Separate gene info (but keep SNP as linking key)
      gene_info_cols <- c("phenotype_id", "variant_id", "GeneSymbol", "EnsemblID", "GeneId")
      dt_for_munge <- processed_dt[, setdiff(names(processed_dt), gene_info_cols), with = FALSE]
      gene_info <- processed_dt[, intersect(names(processed_dt), c("SNP", gene_info_cols)), with = FALSE]

      # --- File names
      input_name <- gsub("\\.tsv$", "", basename(f))
      gene_info_file <- file.path(output_path, paste0("gene_info_", input_name, ".tsv.gz"))
      output_file <- file.path(output_path, paste0("munged_", input_name, ".tsv.gz"))
      merged_output_file <- file.path(output_path, paste0("munged_with_genes_", input_name, ".tsv.gz"))

      # --- Save gene info
      fwrite(gene_info, gene_info_file, sep = "\t")

      # --- Save temp file for MungeSumstats
      temp_file <- tempfile(fileext = ".tsv")
      fwrite(dt_for_munge, temp_file, sep = "\t")

      # --- Run MungeSumstats
      MungeSumstats::format_sumstats(
        path = temp_file,
        ref_genome = "GRCh38",
        convert_ref_genome = "GRCh37",
        check_dups = FALSE,
        compute_z = 'BETA',
        flip_frq_as_biallelic = TRUE,
        bi_allelic_filter = FALSE,
        drop_indels = TRUE,
        save_path = output_file,
        return_data = FALSE,
        log_folder = logs_path,
        log_folder_ind = TRUE,
        log_mungesumstats_msgs = TRUE,
        force_new = TRUE,
        sort_coordinates = TRUE
      )

      # --- Merge munged data with gene info if requested
      if (keep_gene_info_merged && file.exists(output_file)) {
        merge_munged_with_gene_info(output_file, gene_info_file, merged_output_file)
      }

      unlink(temp_file)
      cat("✓ Done:", basename(f), "\n")
      cat("  -> Munged:", basename(output_file), "\n")
      cat("  -> Gene info:", basename(gene_info_file), "\n")
      if (keep_gene_info_merged) {
        cat("  -> Merged:", basename(merged_output_file), "\n")
      }
      cat("\n")

    }, error = function(e) {
      cat("✗ Error with", basename(f), ":", e$message, "\n")
    })
  }

  # --- Summary
  cat("\nMunging complete!\n")
  cat("Munged files:", length(list.files(output_path, pattern = "^munged_.*\\.tsv\\.gz$")), "\n")
  cat("Gene info files:", length(list.files(output_path, pattern = "^gene_info_.*\\.tsv\\.gz$")), "\n")
  if (keep_gene_info_merged) {
    cat("Merged files:", length(list.files(output_path, pattern = "^munged_with_genes_.*\\.tsv\\.gz$")), "\n")
  }
  cat("Logs in:", logs_path, "\n")
  cat("Output dir:", output_path, "\n")
}

# =============================================
# Helper function to check gene info preservation
# =============================================
check_gene_info_preservation <- function(eqtl_path) {
  output_path <- file.path(eqtl_path, "munged")
  
  # Check a sample merged file
  merged_files <- list.files(output_path, pattern = "^munged_with_genes_.*\\.tsv\\.gz$", full.names = TRUE)
  
  if (length(merged_files) > 0) {
    sample_file <- merged_files[1]
    dt <- fread(sample_file, nrows = 10)
    
    cat("Sample of merged file:", basename(sample_file), "\n")
    cat("Columns:", paste(names(dt), collapse = ", "), "\n")
    cat("EnsemblID preserved:", "EnsemblID" %in% names(dt), "\n")
    
    if ("EnsemblID" %in% names(dt)) {
      cat("Sample EnsemblID values:", paste(head(unique(dt$EnsemblID)), collapse = ", "), "\n")
    }
  } else {
    cat("No merged files found. Check if keep_gene_info_merged = TRUE\n")
  }
}

# =============================================
# Example Run
# =============================================
munge_all_fastqtl(
  eqtl_path = "/home/ubuntu/results_eqtl_2592025/",
  pattern = "*.tsv",
  keep_gene_info_merged = TRUE  # This will create files with gene info included
)

# Check that gene info was preserved
check_gene_info_preservation("/home/ubuntu/results_eqtl_2592025/")
