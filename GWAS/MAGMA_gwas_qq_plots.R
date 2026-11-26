# GWAS MAGMA Analysis Script for Specific Diseases
# Load required libraries
library(readr)
library(dplyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(gridExtra)
library(scales)
library(viridis)  # for colorblind-friendly palette

# Define diseases and their labels
diseases <- c("ASD","ADHD", "DD","TS","EPI", "GEPI", "CEPI", "SCZ","MDD","AN","OCD",
              "PANIC", "BIP","PTSD", "PAU","SUD","ALS","FTLD", "MS","PD","ALZ",
              "IQ", "EA", "AOW", "Height")

diseases_labels <- c("Autism","ADHD","Developmental disorders","Tourette's syndrome", "Epilepsy", 
                     "Generalized epilepsy", "Childhood epilepsy","Schizophrenia", "Major depression","Eating disorder", 
                     "OCD","Panic disorder", "Bipolar disorder", "PTSD","Alcohol use","Substance use", 
                     "Amyotrophic lateral sclerosis","Frontotemporal lobar degeneration", "Multiple sclerosis", 
                     "Parkinson's disease", "Alzheimer's disease",
                     "IQ","Educational attainment","Age at onset of walking", "Height")

# Create mapping between disease codes and labels
disease_label_map <- setNames(diseases_labels, diseases)

# Function to load all MAGMA data for a disease
load_disease_magma_data <- function(disease) {
  file_path <- paste0("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Enrichment_analysis/gwas_magma/", disease, "/magma.genes.out")
  
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  tryCatch({
    magma_df <- read_table(file_path)
    magma_df$disease <- disease
    magma_df$disease_label <- disease_label_map[disease]
    return(magma_df)
  }, error = function(e) {
    warning(paste("Error reading file for disease", disease, ":", e$message))
    return(NULL)
  })
}

# Function to analyze all diseases
analyze_all_diseases <- function() {
  cat("Analyzing diseases:", paste(diseases, collapse = ", "), "\n")
  
  all_data <- map_dfr(diseases, load_disease_magma_data)
  
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No data could be loaded")
  }
  
  gene_counts <- all_data %>%
    group_by(disease, disease_label) %>%
    summarise(
      total_genes = n(),
      genes_p_lt_0.05 = sum(P < 0.05, na.rm = TRUE),
      genes_p_lt_0.01 = sum(P < 0.01, na.rm = TRUE),
      genes_p_lt_0.005 = sum(P < 0.005, na.rm = TRUE),
      genes_p_lt_0.001 = sum(P < 0.001, na.rm = TRUE),
      median_p = median(P, na.rm = TRUE),
      mean_p = mean(P, na.rm = TRUE),
      min_p = min(P, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(genes_p_lt_0.05))
  
  disease_gene <- setNames(lapply(diseases, function(disease) {
    tryCatch({
      toupper(load_disease_magma_data(disease))
    }, error = function(e) {
      warning(paste("Could not load genes for", disease, ":", e$message))
      return(character(0))
    })
  }), diseases)
  
  return(list(
    raw_data = all_data,
    summary = gene_counts,
    diseases = diseases,
    disease_genes = disease_gene
  ))
}

# Gene count visualization
plot_gene_counts <- function(analysis_results) {
  gene_counts <- analysis_results$summary
  print(gene_counts$disease)
  print(gene_counts$genes_p_lt_0.05)
  plot_data <- gene_counts %>%
    select(disease, disease_label, total_genes, genes_p_lt_0.05, genes_p_lt_0.01, genes_p_lt_0.005, genes_p_lt_0.001) %>%
    pivot_longer(cols = starts_with("genes") | starts_with("total"), 
                 names_to = "category", values_to = "count") %>%
    mutate(
      category = case_when(
        category == "total_genes" ~ "Total genes",
        category == "genes_p_lt_0.05" ~ "p < 0.05",
        category == "genes_p_lt_0.01" ~ "p < 0.01",
        category == "genes_p_lt_0.005" ~ "p < 0.005",
        category == "genes_p_lt_0.001" ~ "p < 0.001"
      ),
      category = factor(category, levels = c("Total genes", "p < 0.05", "p < 0.01", "p < 0.005", "p < 0.001"))
    )
  
  p1 <- ggplot(plot_data, aes(x = reorder(disease_label, count), y = count, fill = category)) +
    geom_col(position = position_dodge(width = 0.8)) +
    # Only add labels for "p < 0.05"
    geom_text(
      data = filter(plot_data, category == "p < 0.05"),
      aes(label = count),
      position = position_dodge(width = 0.8),
      hjust = -0.1,
      size = 3
    ) +
    coord_flip() +
    scale_fill_viridis_d(option = "D") +
    labs(
      title = "Gene Counts by Disease and Significance Threshold",
      x = "Disease",
      y = "Number of genes by threshold: 0.05",
      fill = "Significance"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.y = element_text(size = 10)
    )
  
  return(p1)
}

# P-value distribution
plot_pvalue_distributions <- function(analysis_results, ncol = 5) {
  all_data <- analysis_results$raw_data
  
  p1 <- ggplot(all_data, aes(x = P)) +
    geom_histogram(bins = 50, fill = viridis(1), alpha = 0.8, color = "black") +
    facet_wrap(~disease_label, scales = "free_y", ncol = ncol) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "P-value Distributions by Disease",
      x = "P-value",
      y = "Gene Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 9, face = "bold")
    )
  
  p2 <- all_data %>%
    filter(P > 0) %>%
    ggplot(aes(x = -log10(P))) +
    geom_histogram(bins = 50, fill = viridis(1, option = "C"), alpha = 0.8, color = "black") +
    facet_wrap(~disease_label, scales = "free_y", ncol = ncol) +
    labs(
      title = "P-value Distributions (Log Scale) by Disease",
      x = expression(-log[10](P)),
      y = "Gene Count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 9, face = "bold")
    )
  
  return(list(linear = p1, log_scale = p2))
}

# QQ plots
plot_qq_plots <- function(analysis_results, ncol = 5) {
  all_data <- analysis_results$raw_data
  
  qq_data <- all_data %>%
    filter(P > 0, !is.na(P)) %>%
    group_by(disease, disease_label) %>%
    arrange(P) %>%
    mutate(
      expected = -log10(ppoints(n())),
      observed = -log10(P)
    ) %>%
    ungroup()
  
  p1 <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6, size = 0.8, color = viridis(1)) +
    geom_abline(slope = 1, intercept = 0, color = "gray30", linetype = "dashed") +
    facet_wrap(~disease_label, scales = "free", ncol = ncol) +
    labs(
      title = "Q-Q Plots by Disease",
      x = "Expected -log10(P)",
      y = "Observed -log10(P)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      strip.text = element_text(size = 9, face = "bold")
    )
  
  return(p1)
}

# Comparison plot
plot_gene_count_comparison <- function(analysis_results) {
  gene_counts <- analysis_results$summary
  
  p1 <- ggplot(gene_counts, aes(x = reorder(disease_label, genes_p_lt_0.05), y = genes_p_lt_0.05)) +
    geom_col(fill = viridis(1), alpha = 0.8) +
    geom_text(
      aes(label = genes_p_lt_0.05),
      hjust = -0.1,
      size = 3
    ) +
    coord_flip() +
    labs(
      title = "Number of Significant Genes (p < 0.05) by Disease",
      x = "Disease",
      y = "Number of Genes",
      subtitle = "Filtered at p < 0.05"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 10)
    ) +
    expand_limits(y = max(gene_counts$genes_p_lt_0.05) * 1.1)  # Add some padding for labels
  
  return(p1)
}

# Print summary table
print_summary_table <- function(analysis_results) {
  summary_table <- analysis_results$summary
  
  cat("\n=== GWAS MAGMA Analysis Summary ===\n")
  cat("Total diseases analyzed:", nrow(summary_table), "\n\n")
  print(summary_table, n = Inf)
  
  cat("\n=== Top 10 Diseases by Significant Genes (p < 0.005) ===\n")
  top10 <- summary_table %>%
    arrange(desc(genes_p_lt_0.005)) %>%
    slice_head(n = 10)
  print(top10)
  
  cat("\n=== Disease Gene List Sizes (p < 0.005, top 500) ===\n")
  gene_sizes <- sapply(analysis_results$disease_genes, length)
  df <- data.frame(
    disease = names(gene_sizes),
    disease_label = disease_label_map[names(gene_sizes)],
    gene_count = gene_sizes
  ) %>% arrange(desc(gene_count))
  print(df)
}

# Main analysis
main_analysis <- function() {
  cat("Starting GWAS MAGMA analysis for specified diseases...\n")
  
  results <- analyze_all_diseases()
  print_summary_table(results)
  
  cat("\nCreating visualizations...\n")
  
  p1 <- plot_gene_counts(results)
  p2 <- plot_gene_count_comparison(results)
  pval_plots <- plot_pvalue_distributions(results)
  qq <- plot_qq_plots(results)
  
  print(p1)
  print(p2)
  print(pval_plots$linear)
  print(pval_plots$log_scale)
  print(qq)
  
  ggsave("gene_counts_by_disease.pdf", p1, width = 7, height = 4)
  ggsave("gene_comparison_p005.pdf", p2, width = 7, height = 4)
  ggsave("pvalue_distributions.pdf", pval_plots$linear, width = 18, height = 12)
  ggsave("pvalue_distributions_log.pdf", pval_plots$log_scale, width = 18, height = 12)
  ggsave("qq_plots.png", qq, width = 18, height = 12)
  
  cat("\nAll plots saved to current directory.\n")
  
  return(results)
}

# Run analysis
results <- main_analysis()
print_summary_table(results)
plot_gene_counts(results)