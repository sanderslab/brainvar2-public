library(tidyverse)
library(coloc)
library(data.table)

run_coloc <- function(input_dir, N_gwas, N_eqtl, output_file) {
  
  all_results <- list()
  
  gwas_files <- dir(input_dir, pattern = ".*gwas.txt", full.names = TRUE)
  eqtl_files <- dir(input_dir, pattern = ".*eqtl.txt", full.names = TRUE)
  
  cat("Found", length(gwas_files), "GWAS files and", length(eqtl_files), "eQTL files\n")
  
  for (i in seq_along(gwas_files)) {
    cat("\n=== Processing file", i, "of", length(gwas_files), "===\n")
    
    tryCatch({
      # Read data
      gwas <- read_tsv(gwas_files[i], show_col_types = FALSE)
      eqtl <- read_tsv(eqtl_files[i], show_col_types = FALSE)
      
      # Standardize column types
      gwas <- gwas %>%
        mutate(
          SNP = as.character(SNP),
          CHR = as.character(CHR),
          BP = as.numeric(BP),
          A1 = as.character(A1),
          A2 = as.character(A2)
        )
      
      eqtl <- eqtl %>%
        mutate(
          SNP = as.character(SNP),
          CHR = as.character(CHR), 
          BP = as.numeric(BP),
          A1 = as.character(A1),
          A2 = as.character(A2)
        )
      
      if (nrow(gwas) == 0 || nrow(eqtl) == 0) next
      
      # Handle BETA: convert OR to log(OR) if needed
      # This is the standard and statistically correct approach for colocalization
      if (!"BETA" %in% names(gwas) && "OR" %in% names(gwas)) {
         cat("Converting OR to BETA using log transformation (standard approach)\n")
         gwas <- gwas %>% mutate(BETA = log(OR))
        }

      # Handle SE: estimate if missing
      if (!"SE" %in% names(gwas)) {
       cat("Estimating SE from available data\n")

       # Method 1: From P-value and BETA
      if ("P" %in% names(gwas) && "BETA" %in% names(gwas)) {
          gwas <- gwas %>%
          mutate(
           Z = sign(BETA) * abs(qnorm(P / 2, lower.tail = FALSE)),
           SE = abs(BETA / Z)
         )

        # Method 3: From OR and P-value
        } else if ("OR" %in% names(gwas) && "P" %in% names(gwas)) {
         gwas <- gwas %>%
         mutate(
         BETA = log(OR),
         Z = sign(BETA) * abs(qnorm(P / 2, lower.tail = FALSE)),
         SE = abs(BETA / Z)
           )

        # Method 4: From confidence intervals
         } else if ("OR_L95" %in% names(gwas) && "OR_U95" %in% names(gwas)) {
         gwas <- gwas %>%
         mutate(
         BETA = log(OR),
         SE = (log(OR_U95) - log(OR_L95)) / 3.92
         )
         }
    # Method 2: From N and Z
        else if ("N" %in% names(gwas) && "Z" %in% names(gwas)) {
         gwas <- gwas %>%
         mutate(
         N = as.numeric(N),
         SE = 1 / sqrt(N),
         BETA = Z * SE
         )
        }
  }



















      
      # Rename columns
      gwas <- gwas %>% rename(BETA_gwas = BETA, SE_gwas = SE)
      eqtl <- eqtl %>% rename(BETA_eqtl = BETA, SE_eqtl = SE)
      
      # Rename P columns safely
      if ("P" %in% names(gwas)) gwas <- gwas %>% rename(P_gwas = P)
      if ("pvalue" %in% names(gwas)) gwas <- gwas %>% rename(P_gwas = pvalue)
      if ("PVAL" %in% names(gwas)) gwas <- gwas %>% rename(P_gwas = PVAL)
      
      if ("P" %in% names(eqtl)) eqtl <- eqtl %>% rename(P_eqtl = P)
      if ("pvalue" %in% names(eqtl)) eqtl <- eqtl %>% rename(P_eqtl = pvalue)
      if ("PVAL" %in% names(eqtl)) eqtl <- eqtl %>% rename(P_eqtl = PVAL)
      
      # Create variant IDs
      gwas <- gwas %>% mutate(variant_id = paste(SNP, A1, A2, sep = "_"))
      eqtl <- eqtl %>% mutate(variant_id = paste(SNP, A1, A2, sep = "_"))
      
      # Handle FRQ: use eQTL FRQ if GWAS FRQ is missing
      if (!"FRQ" %in% names(gwas)) {
        cat("Using eQTL FRQ for GWAS (FRQ not in GWAS file)\n")
        gwas <- gwas %>%
          left_join(eqtl %>% select(variant_id, FRQ), by = "variant_id")
      }
      
      gwas <- gwas %>% mutate(FRQ_gwas = FRQ)
      
      # Use eQTL FRQ
      if ("FRQ" %in% names(eqtl)) {
        eqtl <- eqtl %>% mutate(FRQ_eqtl = FRQ)
      } else {
        stop("eQTL file missing FRQ column")
      }
      
      # Join on variant
      overlap <- inner_join(
        gwas, eqtl,
        by = c("variant_id", "SNP", "CHR", "BP", "A1", "A2"),
        suffix = c("_gwas", "_eqtl")
      )
      
      if (nrow(overlap) == 0) next
      
      # Get lead SNPs
      lead_gwas <- overlap %>%
        arrange(P_gwas) %>%
        slice(1) %>%
        select(variant_id, BETA_gwas, P_gwas, BETA_eqtl, P_eqtl)
      
      lead_eqtl <- overlap %>%
        arrange(P_eqtl) %>%
        slice(1) %>%
        select(variant_id, BETA_gwas, P_gwas, BETA_eqtl, P_eqtl)
      
      # Prepare GWAS data
      gwas_clean <- overlap %>%
        transmute(
          beta = BETA_gwas,
          varbeta = SE_gwas^2,
          snp = variant_id,
          pvalues = P_gwas,
          position = BP,
          MAF = pmin(FRQ_gwas, 1 - FRQ_gwas)
        ) %>% filter(complete.cases(.))
      
      # Prepare eQTL data
      eqtl_clean <- overlap %>%
        transmute(
          beta = BETA_eqtl,
          varbeta = SE_eqtl^2,
          snp = variant_id,
          pvalues = P_eqtl,
          position = BP,
          MAF = pmin(FRQ_eqtl, 1 - FRQ_eqtl)
        ) %>% filter(complete.cases(.))
      
      if (nrow(gwas_clean) == 0 || nrow(eqtl_clean) == 0) next
      
      # Prepare inputs for coloc
      gwas_input <- as.list(gwas_clean)
      gwas_input$N <- N_gwas
      gwas_input$type <- 'quant'
      
      eqtl_input <- as.list(eqtl_clean)
      eqtl_input$N <- N_eqtl
      eqtl_input$type <- 'quant'
      
      # Run coloc
      res <- coloc.abf(dataset1 = gwas_input, dataset2 = eqtl_input)
      
      # Collect gene info
      gene_info <- tibble(
        chr = paste0("chr", unique(gwas$CHR)),
        gene_ensembl = if ("ENSEMBLID" %in% names(eqtl)) unique(eqtl$ENSEMBLID) else NA_character_,
        gene_name = if ("GENESYMBOL" %in% names(eqtl)) unique(eqtl$GENESYMBOL) else NA_character_
      )
      
      # Combine results with lead SNPs
      final_res <- cbind(
        gene_info,
        as.data.frame(t(res$summary))[, c(
          "nsnps", "PP.H0.abf", "PP.H1.abf",
          "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"
        )],
        gwas_lead_variant = lead_gwas$variant_id,
        gwas_lead_beta = lead_gwas$BETA_gwas,
        gwas_lead_pval = lead_gwas$P_gwas,
        gwas_lead_beta_eqtl = lead_gwas$BETA_eqtl,
        gwas_lead_pval_eqtl = lead_gwas$P_eqtl,
        eqtl_lead_variant = lead_eqtl$variant_id,
        eqtl_lead_beta_gwas = lead_eqtl$BETA_gwas,
        eqtl_lead_pval_gwas = lead_eqtl$P_gwas,
        eqtl_lead_beta_eqtl = lead_eqtl$BETA_eqtl,
        eqtl_lead_pval_eqtl = lead_eqtl$P_eqtl
      )
      
      all_results <- append(all_results, list(final_res))
      cat("✓ Processed gene:", gene_info$gene_name, "\n")
      
    }, error = function(e) {
      cat("✗ Error in file", i, ":", e$message, "\n")
    })
  }
  
  if (length(all_results) == 0) {
    cat("No results generated.\n")
    return(NULL)
  }
  
  all_results_df <- bind_rows(all_results)
  
  # Save results
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  write_tsv(all_results_df, output_file)
  
  cat("\n=== COLOC ANALYSIS COMPLETE ===\nResults saved to:", output_file, "\n")
  return(all_results_df)
}

# Define GWAS configurations
gwas_configs <- list(
  AOW = list(name = "AOW", N = 350000),
  ASD = list(name = "ASD", N = 46350),
  SCZ = list(name = "SCZ", N = 320404),
  ALZ = list(name = "ALZ", N = 487511),
  AN = list(name = "AN", N = 72517),
  EA = list(name = "EA", N = 765283),
  EPI = list(name = "EPI", N = 45000),
  IQ = list(name = "IQ", N = 269867),
  MDD = list(name = "MDD", N = 1564598),
  MS = list(name = "MS", N = 41505),
  OCD = list(name = "OCD", N = 9725),
  PD = list(name = "PD", N = 611485),
  TS = list(name = "TS", N = 14307),
  PAU = list(name = "PAU", N = 1079947),
  BIP = list(name = "BIP", N = 413466),
  PANIC = list(name = "PANIC", N = 23000),
  ADHD = list(name = "ADHD", N = 225534),
  DD = list(name = "DD", N = 32952),
  GEPI = list(name = "GEPI", N = 59945),
  CEPI = list(name = "CEPI", N = 43485),
  Height = list(name = "Height", N = 5400000),
  SUD = list(name = "SUD", N = 1025550),
  PTSD = list(name = "PTSD", N = 1222882),
  OCD2025 = list(name = "OCD2025", N = 53234),
  FTLD_TDP = list(name = "FTLD_TDP", N = 7327),
  ALS = list(name = "ALS", N = 138086)
)
gwas_configs <- list(
PTSD = list(name = "PTSD", N = 1222882)
)

# Define eQTL configurations
eqtl_configs <- list(
  bmind_postnatal_ExN = 113,
  bmind_postnatal_InN = 113,
  bmind_postnatal_Glia = 113,
  bmind_postnatal_NNC = 113,
  bmind_prenatal_ExN = 150,
  bmind_prenatal_InN = 150,
  bmind_prenatal_Glia = 150,
  bmind_prenatal_NNC = 150,
  pseudobulk_prenatal_ExN = 47,
  pseudobulk_prenatal_InN = 47,
  pseudobulk_prenatal_Glia = 47,
  pseudobulk_prenatal_NNC = 47,
  pseudobulk_postnatal_ExN = 71,
  pseudobulk_postnatal_InN = 71,
  pseudobulk_postnatal_Glia = 71,
  pseudobulk_postnatal_NNC = 71,
  bulk_postnatal = 150,
  bulk_prenatal = 150
)

# Run all combinations
run_all_coloc <- function(gwas_list = names(gwas_configs), 
                          eqtl_list = names(eqtl_configs),
                          base_dir = "/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/intersections") {
  
  results_summary <- list()
  
  for (gwas_name in gwas_list) {
    for (eqtl_name in eqtl_list) {
      
      cat("\n\n========================================\n")
      cat("Processing:", gwas_name, "x", eqtl_name, "\n")
      cat("========================================\n")
      
      input_dir <- file.path(base_dir, paste0(gwas_name, "_", eqtl_name, "/"))
      output_dir <- file.path(base_dir, paste0(gwas_name, "_", eqtl_name, "_targeted/"))
      output_file <- '/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/'
      output_file <- file.path(output_dir, paste0(gwas_name, "_", eqtl_name, "_targeted.txt"))
      
      if (!dir.exists(input_dir)) {
        cat("Input directory not found, skipping:", input_dir, "\n")
        next
      }
      
      result <- run_coloc(
        input_dir = input_dir,
        N_gwas = gwas_configs[[gwas_name]]$N,
        N_eqtl = eqtl_configs[[eqtl_name]],
        output_file = output_file
      )
      
      if (!is.null(result)) {
        results_summary[[paste0(gwas_name, "_", eqtl_name)]] <- list(
          gwas = gwas_name,
          eqtl = eqtl_name,
          n_genes = nrow(result),
          n_coloc = sum(result$PP.H4.abf > 0.8, na.rm = TRUE),
          mean_pp4 = mean(result$PP.H4.abf, na.rm = TRUE),
          output = output_file
        )
      }
    }
  }
  
  cat("\n\n========================================\n")
  cat("ALL ANALYSES COMPLETE\n")
  cat("========================================\n")
  cat("Total combinations processed:", length(results_summary), "\n")
  
  return(results_summary)
}

# ============================================================================
# TEST RUN: SCZ ONLY
# ============================================================================

run_scz_test <- function(base_dir = "/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/intersections") {
  
  cat("========================================\n")
  cat("COLOCALIZATION ANALYSIS TEST RUN: SCZ\n")
  cat("========================================\n\n")
  
  # SCZ configuration
  gwas_name <- "SCZ"
  N_gwas <- 320404
  
  # Store results summary
  results_summary <- list()
  
  # Run SCZ colocalization for all eQTL datasets
  for (eqtl_name in names(eqtl_configs)) {
    
    cat("\n========================================\n")
    cat("Processing: SCZ x", eqtl_name, "\n")
    cat("========================================\n")
    
    # Define paths
    input_dir <- file.path(base_dir, paste0(gwas_name, "_", eqtl_name, "/"))
    output_dir <- file.path(base_dir, paste0(gwas_name, "_", eqtl_name, "_targeted/"))
    output_file <- file.path(output_dir, paste0(gwas_name, "_", eqtl_name, "_targeted.txt"))
    
    # Check if input directory exists
    if (!dir.exists(input_dir)) {
      cat("⚠ Input directory not found, skipping:", input_dir, "\n")
      next
    }
    
    # Run colocalization
    tryCatch({
      result <- run_coloc(
        input_dir = input_dir,
        N_gwas = N_gwas,
        N_eqtl = eqtl_configs[[eqtl_name]],
        output_file = output_file
      )
      
      if (!is.null(result)) {
        # Store summary
        results_summary[[eqtl_name]] <- list(
          eqtl = eqtl_name,
          n_genes = nrow(result),
          n_coloc = sum(result$PP.H4.abf > 0.8, na.rm = TRUE),
          mean_pp4 = mean(result$PP.H4.abf, na.rm = TRUE),
          output = output_file
        )
        
        cat("✓ Successfully processed:", eqtl_name, "\n")
        cat("  - Genes analyzed:", nrow(result), "\n")
        cat("  - Strong colocalizations (PP.H4 > 0.8):", 
            sum(result$PP.H4.abf > 0.8, na.rm = TRUE), "\n")
      }
      
    }, error = function(e) {
      cat("✗ Error processing", eqtl_name, ":", e$message, "\n")
    })
  }
  
  # Print final summary
  cat("\n\n========================================\n")
  cat("SCZ COLOCALIZATION TEST COMPLETE\n")
  cat("========================================\n\n")
  
  if (length(results_summary) > 0) {
    summary_df <- bind_rows(results_summary)
    
    cat("Summary of results:\n")
    cat("-------------------\n")
    for (i in seq_along(results_summary)) {
      res <- results_summary[[i]]
      cat(sprintf("%-30s: %3d genes, %3d coloc (PP.H4>0.8), mean PP.H4=%.3f\n",
                  res$eqtl, res$n_genes, res$n_coloc, res$mean_pp4))
    }
    
    cat("\n")
    cat("Total eQTL datasets processed:", length(results_summary), "\n")
    cat("Total genes analyzed:", sum(sapply(results_summary, function(x) x$n_genes)), "\n")
    cat("Total strong colocalizations:", sum(sapply(results_summary, function(x) x$n_coloc)), "\n")
    
    # Save summary
    summary_file <- file.path(base_dir, "SCZ_coloc_summary.txt")
    write_tsv(summary_df, summary_file)
    cat("\nSummary saved to:", summary_file, "\n")
    
  } else {
    cat("⚠ No results generated. Please check:\n")
    cat("  1. Input directory exists and contains data\n")
    cat("  2. File naming conventions match expected patterns\n")
    cat("  3. Required columns are present in input files\n")
  }
  
  cat("\n========================================\n")
  
  return(results_summary)
}

# ============================================================================
# USAGE EXAMPLES
# ============================================================================

# Example 1: Run SCZ test only
#scz_results <- run_scz_test()

# Example 2: Run all analyses
 all_results <- run_all_coloc()

# Example 3: Run specific GWAS traits
# results <- run_all_coloc(gwas_list = c("TS", "ASD", "SCZ"))

# Example 4: Run specific eQTL types
# results <- run_all_coloc(eqtl_list = c("bmind_postnatal_ExN", "bmind_prenatal_ExN"))

# Example 5: Run single analysis manually
# result <- run_coloc(
#   input_dir = '/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/intersections/SCZ_bmind_postnatal_ExN/',
#   N_eqtl = 113,
#   N_gwas = 320404,
#   output_file = '/mnt/s3/eqtl/eQTL/colocalisation_eqtl_subset/intersections/SCZ_bmind_postnatal_ExN_targeted/SCZ_bmind_postnatal_ExN_targeted.txt'
# )
