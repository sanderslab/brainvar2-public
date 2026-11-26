# ============================================================================
# RNA GENE ENRICHMENT PIPELINE
# ============================================================================

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidyr)
library(purrr)
library(MAGMA.Celltyping)
library(MungeSumstats)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.1000genomes.hs37d5)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(fgsea)
library(openxlsx)

set.seed(123)
setwd("/Users/nrohani/Documents/Sanders_lab/BrainVar/Gene-enrichment/Enrichment_analysis")

# ============================================================================
# SHARED HELPERS
# ============================================================================

save_enrichment_results <- function(enrichment_scores_matrix,
                                    p_values_matrix,
                                    adj_p_values_matrix,
                                    fgsea_results_list,
                                    output_csv_dir,
                                    label) {
  if (!dir.exists(output_csv_dir)) dir.create(output_csv_dir, recursive = TRUE)
  
  write.csv(enrichment_scores_matrix, file.path(output_csv_dir, paste0(label, "_NES.csv")))
  write.csv(p_values_matrix,         file.path(output_csv_dir, paste0(label, "_pval.csv")))
  write.csv(adj_p_values_matrix,     file.path(output_csv_dir, paste0(label, "_adjpval.csv")))
  
  leading_edge_rows <- list()
  for (group_name in names(fgsea_results_list)) {
    res_df <- fgsea_results_list[[group_name]]
    if (is.null(res_df) || nrow(res_df) == 0) next
    
    for (i in seq_len(nrow(res_df))) {
      genes <- res_df$leadingEdge[[i]]
      leading_edge_rows[[length(leading_edge_rows) + 1]] <- data.frame(
        Group            = group_name,
        Disease          = res_df$pathway[i],
        NES              = res_df$NES[i],
        pval             = res_df$pval[i],
        padj             = res_df$padj[i],
        LeadingEdgeGenes = if (length(genes) > 0) paste(genes, collapse = ", ") else NA_character_,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(leading_edge_rows) > 0) {
    leading_edge_df <- bind_rows(leading_edge_rows)
    write.csv(leading_edge_df,
              file.path(output_csv_dir, paste0(label, "_leading_edge_genes.csv")),
              row.names = FALSE)
  }
  
  cat("  Saved CSVs to:", output_csv_dir, "with prefix:", label, "\n")
}

flatten_fgsea_results <- function(fgsea_results_list,
                                  cell_type = NULL,
                                  stage     = NULL,
                                  modality  = "") {
  rows <- list()
  for (group_name in names(fgsea_results_list)) {
    res_df <- fgsea_results_list[[group_name]]
    if (is.null(res_df) || nrow(res_df) == 0) next
    
    for (i in seq_len(nrow(res_df))) {
      genes <- res_df$leadingEdge[[i]]
      rows[[length(rows) + 1]] <- data.frame(
        CellSubtype = if (!is.null(cell_type)) cell_type else group_name,
        Stage       = if (!is.null(stage)) stage else group_name,
        Disease     = res_df$pathway[i],
        NES         = res_df$NES[i],
        pval        = res_df$pval[i],
        padj        = res_df$padj[i],
        LeadingEdge = if (length(genes) > 0) paste(genes, collapse = ", ") else NA_character_,
        Modality    = modality,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

append_sheet <- function(wb, tidy_df, sheet_name) {
  sheet_name <- substr(sheet_name, 1, 31)
  
  hs <- createStyle(
    fontColour     = "#FFFFFF",
    fgFill         = "#2E4057",
    halign         = "CENTER",
    textDecoration = "Bold",
    border         = "Bottom"
  )
  
  addWorksheet(wb, sheet_name)
  writeDataTable(wb, sheet = sheet_name, x = tidy_df,
                 tableStyle = "TableStyleMedium2", withFilter = TRUE)
  addStyle(wb, sheet_name, style = hs,
           rows = 1, cols = seq_len(ncol(tidy_df)), gridExpand = TRUE)
  
  setColWidths(wb, sheet_name, cols = seq_len(ncol(tidy_df)), widths = "auto")
  le_col <- which(colnames(tidy_df) == "LeadingEdge")
  if (length(le_col)) setColWidths(wb, sheet_name, cols = le_col, widths = 60)
}

# ============================================================================
# DISEASE ORDERING / GROUPING (ATAC-STYLE)
# ============================================================================

DISEASE_ORDER <- c(
  "NDD_rare", "DD", "ASD_rare", "ASD", "ADHD",
  "CEPI", "GEPI", "TS", "AN", "PAU", "SUD",
  "SCZ", "BIP", "PANIC", "OCD", "PTSD", "MDD",
  "MS", "ALS", "FTLD", "PD", "ALZ",
  "Height", "AOW", "IQ", "EA"
)

DISEASE_LABELS <- c(
  "Developmental disorders rare", "Developmental disorders",
  "Autism rare", "Autism", "ADHD",
  "Childhood epilepsy", "Generalized epilepsy", "Tourette's syndrome",
  "Eating disorder", "Alcohol use", "Substance use",
  "Schizophrenia", "Bipolar disorder", "Panic disorder",
  "OCD", "PTSD", "Major depression",
  "Multiple sclerosis", "Amyotrophic lateral sclerosis",
  "Frontotemporal lobar degeneration", "Parkinson's disease",
  "Alzheimer's disease",
  "Height", "Age at onset of walking", "IQ", "Educational attainment"
)

DISEASE_CATEGORIES <- c(
  "NDD_rare" = "Neurodevelopmental", "DD"    = "Neurodevelopmental",
  "ASD_rare" = "Neurodevelopmental", "ASD"   = "Neurodevelopmental",
  "ADHD"     = "Neurodevelopmental", "CEPI"  = "Neurodevelopmental",
  "GEPI"     = "Neurodevelopmental", "TS"    = "Neurodevelopmental",
  "AN"       = "Neuropsychiatric",   "PAU"   = "Neuropsychiatric",
  "SUD"      = "Neuropsychiatric",   "SCZ"   = "Neuropsychiatric",
  "BIP"      = "Neuropsychiatric",   "PANIC" = "Neuropsychiatric",
  "OCD"      = "Neuropsychiatric",   "PTSD"  = "Neuropsychiatric",
  "MDD"      = "Neuropsychiatric",
  "MS"       = "Neurodegenerative",  "ALS"   = "Neurodegenerative",
  "FTLD"     = "Neurodegenerative",  "PD"    = "Neurodegenerative",
  "ALZ"      = "Neurodegenerative",
  "Height"   = "Cognitive / other",  "AOW"   = "Cognitive / other",
  "IQ"       = "Cognitive / other",  "EA"    = "Cognitive / other"
)

CATEGORY_LEVELS <- c("Neurodevelopmental", "Neuropsychiatric", "Neurodegenerative", "Cognitive / other")

reorder_disease_rows <- function(mat) {
  label_map     <- setNames(DISEASE_LABELS, DISEASE_ORDER)
  existing_rows <- intersect(DISEASE_ORDER, rownames(mat))
  mat           <- mat[existing_rows, , drop = FALSE]
  rownames(mat) <- label_map[existing_rows]
  row_split     <- factor(DISEASE_CATEGORIES[existing_rows], levels = CATEGORY_LEVELS)
  
  list(
    mat           = mat,
    label_map     = label_map,
    existing_rows = existing_rows,
    row_split     = row_split
  )
}

# ============================================================================
# DATA PREP HELPERS
# ============================================================================

load_marker_genes <- function(directory) {
  deg_files <- list.files(directory,
                          pattern    = "^markers_.*_prenatal_vs_postnatal\\.csv$",
                          full.names = TRUE)
  if (length(deg_files) == 0) stop("No marker files found in: ", directory)
  
  term2gene_list <- list()
  for (file in deg_files) {
    cell_type <- str_extract(basename(file), "(?<=markers_).*(?=_prenatal_vs_postnatal)")
    deg_df    <- read_csv(file, show_col_types = FALSE)
    colnames(deg_df)[1] <- "Gene"
    
    if ("Gene" %in% colnames(deg_df)) {
      filtered <- deg_df %>%
        mutate(Stage = case_when(
          avg_log2FC > 0 & p_val_adj < 1000 ~ "Prenatal",
          avg_log2FC < 0 & p_val_adj < 1000 ~ "Postnatal",
          TRUE ~ NA_character_
        )) %>%
        mutate(term = paste(cell_type, Stage, sep = "_"))
      term2gene_list[[cell_type]] <- filtered
    }
  }
  bind_rows(term2gene_list)
}

prepare_deg_data <- function(deg_data, pct_threshold = 0.1) {
  deg_data$gene <- toupper(deg_data$gene)
  deg_data$gene <- sub("\\..*", "", deg_data$gene)
  if (all(c("pct.1", "pct.2") %in% colnames(deg_data))) {
    deg_data <- deg_data[(deg_data$pct.1 > pct_threshold | deg_data$pct.2 > pct_threshold), ]
  }
  deg_data <- deg_data[!grepl("^ENSG", deg_data$gene), ]
  deg_data
}

select_per_cell <- function(deg_data, cell_types) {
  result_list <- list()
  for (cell_type in cell_types) {
    result_list[[cell_type]] <- deg_data[deg_data$cluster == cell_type, ]
  }
  result_list
}

load_disease_genes <- function(disease, top_n = 500) {
  magma_df <- read_table(paste0("gwas_magma/", disease, "/magma.genes.out"))
  magma_df <- magma_df[magma_df$P < 0.05, ]
  magma_df %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(P) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    slice_head(n = top_n) %>%
    pull(SYMBOL) %>%
    toupper()
}

create_ranked_genes_from_markers <- function(marker_subset) {
  if (!"gene" %in% colnames(marker_subset)) colnames(marker_subset)[1] <- "gene"
  marker_subset$gene <- toupper(marker_subset$gene)
  rankings <- sign(marker_subset$avg_log2FC) * (-log10(marker_subset$p_val_adj))
  names(rankings) <- marker_subset$gene
  rankings <- sort(rankings, decreasing = TRUE)
  max_r <- max(rankings[is.finite(rankings)])
  min_r <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_r, max_r * 10)
  rankings <- replace(rankings, rankings < min_r, min_r * 10)
  sort(rankings, decreasing = TRUE)
}

select_significant_degs <- function(deg_df) {
  deg_df   <- data.frame(deg_df)
  deg_df   <- deg_df[, -1]
  deg_df   <- deg_df[deg_df$p_val_adj < 0.05, ]
  rankings <- sign(deg_df$avg_log2FC) * (-log10(deg_df$p_val_adj))
  names(rankings) <- deg_df$gene
  rankings <- sort(rankings, decreasing = TRUE)
  max_r <- max(rankings[is.finite(rankings)])
  min_r <- min(rankings[is.finite(rankings)])
  rankings <- replace(rankings, rankings > max_r, max_r * 10)
  rankings <- replace(rankings, rankings < min_r, min_r * 10)
  sort(rankings, decreasing = TRUE)
}

# ============================================================================
# MATRIX INTERLEAVING HELPERS (ATAC-STYLE)
# ============================================================================

interleave_pre_post_matrices <- function(mat_pre, mat_post, cell_types) {
  col_labels <- as.vector(rbind(paste0(cell_types, "\nPrenatal"),
                                paste0(cell_types, "\nPostnatal")))
  combined <- matrix(NA, nrow = nrow(mat_pre), ncol = length(col_labels),
                     dimnames = list(rownames(mat_pre), col_labels))
  for (ct in cell_types) {
    if (ct %in% colnames(mat_pre))  combined[, paste0(ct, "\nPrenatal")]  <- mat_pre[, ct]
    if (ct %in% colnames(mat_post)) combined[, paste0(ct, "\nPostnatal")] <- mat_post[, ct]
  }
  combined
}

interleave_pre_post_vector <- function(vec_pre, vec_post, cell_types) {
  out <- as.vector(rbind(vec_pre[cell_types], vec_post[cell_types]))
  names(out) <- as.vector(rbind(paste0(cell_types, "\nPrenatal"),
                                paste0(cell_types, "\nPostnatal")))
  out
}

# ============================================================================
# FGSEA RUNNERS
# ============================================================================

run_fgsea <- function(cell_types, disease_gene, ranked_deg_genes_list, stage,
                      label = "results", output_csv_dir = "enrichment_csv_results") {
  diseases <- names(disease_gene)
  
  enrichment_scores_matrix <- matrix(NA, nrow = length(diseases), ncol = length(cell_types),
                                     dimnames = list(diseases, cell_types))
  p_values_matrix          <- matrix(NA, nrow = length(diseases), ncol = length(cell_types),
                                     dimnames = list(diseases, cell_types))
  adj_p_values_matrix      <- matrix(NA, nrow = length(diseases), ncol = length(cell_types),
                                     dimnames = list(diseases, cell_types))
  
  fgsea_results_list <- list()
  tsv_output_dir <- "fgsea_enriched_genes"
  if (!dir.exists(tsv_output_dir)) dir.create(tsv_output_dir)
  
  for (cell_type in cell_types) {
    deg_genes <- ranked_deg_genes_list[[cell_type]]
    
    enrichment_result <- fgseaMultilevel(
      pathways  = disease_gene,
      stats     = deg_genes,
      minSize   = 5,
      maxSize   = 500,
      scoreType = 'pos',
      eps       = 0.0
    )
    
    cat("\n=== Cell Type:", cell_type, "===\n")
    filtered <- enrichment_result[enrichment_result$NES > 0, ]
    
    if (nrow(filtered) > 0) {
      enriched_list <- purrr::map2_dfr(
        filtered$pathway, filtered$leadingEdge,
        ~data.frame(Disease = .x, LeadingEdgeGene = .y, stringsAsFactors = FALSE)
      )
      out_file <- file.path(tsv_output_dir,
                            paste0("seperate_enriched_genes_", stage, "_",
                                   gsub("[/\\:*?\"<>|]", "_", cell_type), ".tsv"))
      readr::write_tsv(enriched_list, out_file)
      cat("  Saved enriched genes to:", out_file, "\n")
    }
    
    for (i in seq_len(nrow(filtered))) {
      cat("Pathway:", filtered$pathway[i], "\n")
      cat("Leading Edge Genes:\n")
      cat(paste(filtered$leadingEdge[[i]], collapse = ", "), "\n\n")
    }
    
    fgsea_results_list[[cell_type]] <- data.frame(enrichment_result)
    enrichment_result_df <- data.frame(enrichment_result)
    rownames(enrichment_result_df) <- enrichment_result_df$pathway
    
    for (disease in diseases) {
      if (disease %in% rownames(enrichment_result_df)) {
        enrichment_scores_matrix[disease, cell_type] <- enrichment_result_df[disease, ]$NES[1]
        p_values_matrix[disease, cell_type]          <- enrichment_result_df[disease, ]$pval[1]
        adj_p_values_matrix[disease, cell_type]      <- enrichment_result_df[disease, ]$padj[1]
      }
    }
  }
  
  save_enrichment_results(
    enrichment_scores_matrix, p_values_matrix, adj_p_values_matrix,
    fgsea_results_list, output_csv_dir, paste0(label, "_", stage)
  )
  
  tidy_df <- flatten_fgsea_results(
    fgsea_results_list,
    stage    = stage,
    modality = label
  )
  
  list(
    enrichment_scores_matrix = enrichment_scores_matrix,
    p_values_matrix          = p_values_matrix,
    adj_p_values_matrix      = adj_p_values_matrix,
    fgsea_results_list       = fgsea_results_list,
    tidy_df                  = tidy_df
  )
}

run_fgsea_prenatal_postnatal <- function(marker_data, disease_gene, cell_type,
                                         label = "prenatal_postnatal",
                                         output_csv_dir = "enrichment_csv_results") {
  diseases <- names(disease_gene)
  stages   <- c("Prenatal", "Postnatal")
  
  enrichment_scores_matrix <- matrix(NA, nrow = length(diseases), ncol = length(stages),
                                     dimnames = list(diseases, stages))
  p_values_matrix          <- matrix(NA, nrow = length(diseases), ncol = length(stages),
                                     dimnames = list(diseases, stages))
  adj_p_values_matrix      <- matrix(NA, nrow = length(diseases), ncol = length(stages),
                                     dimnames = list(diseases, stages))
  
  fgsea_results_list <- list()
  tsv_output_dir <- "fgsea_enriched_genes"
  if (!dir.exists(tsv_output_dir)) dir.create(tsv_output_dir)
  
  for (stage in stages) {
    stage_data <- marker_data %>% filter(cell_type == !!cell_type, Stage == stage)
    
    if (nrow(stage_data) < 5) {
      cat("  Skipping", stage, "- insufficient genes (", nrow(stage_data), ")\n")
      next
    }
    
    cat("\n  Analyzing", stage, "with", nrow(stage_data), "genes\n")
    stage_data <- prepare_deg_data(stage_data, 0.2)
    rankings   <- create_ranked_genes_from_markers(stage_data)
    
    enrichment_result <- fgseaMultilevel(
      pathways  = disease_gene,
      stats     = rankings,
      minSize   = 5,
      maxSize   = 500,
      scoreType = 'pos',
      eps       = 0.0
    )
    
    filtered <- enrichment_result[enrichment_result$NES > 0, ]
    
    if (nrow(filtered) > 0) {
      enriched_list <- purrr::map2_dfr(
        filtered$pathway, filtered$leadingEdge,
        ~data.frame(Disease = .x, LeadingEdgeGene = .y, stringsAsFactors = FALSE)
      )
      out_file <- file.path(tsv_output_dir,
                            paste0("run_fgsea_prenatal_postnatal_enriched_genes_",
                                   cell_type, "_", stage, ".tsv"))
      readr::write_tsv(enriched_list, out_file)
      cat("  Saved enriched genes to:", out_file, "\n")
    }
    
    if (nrow(filtered) > 0) {
      for (i in seq_len(nrow(filtered))) {
        cat("    Pathway:", filtered$pathway[i], "NES:", round(filtered$NES[i], 3), "\n")
      }
    }
    
    fgsea_results_list[[stage]] <- data.frame(enrichment_result)
    enrichment_result_df <- data.frame(enrichment_result)
    rownames(enrichment_result_df) <- enrichment_result_df$pathway
    
    for (disease in diseases) {
      if (disease %in% rownames(enrichment_result_df)) {
        enrichment_scores_matrix[disease, stage] <- enrichment_result_df[disease, ]$NES[1]
        p_values_matrix[disease, stage]          <- enrichment_result_df[disease, ]$pval[1]
        adj_p_values_matrix[disease, stage]      <- enrichment_result_df[disease, ]$padj[1]
      }
    }
  }
  
  safe_cell_type <- gsub("[/\\:*?\"<>| ]", "_", cell_type)
  save_enrichment_results(
    enrichment_scores_matrix, p_values_matrix, adj_p_values_matrix,
    fgsea_results_list, output_csv_dir, paste0(label, "_", safe_cell_type)
  )
  
  tidy_df <- flatten_fgsea_results(
    fgsea_results_list,
    cell_type = cell_type,
    modality  = label
  )
  
  list(
    enrichment_scores_matrix = enrichment_scores_matrix,
    p_values_matrix          = p_values_matrix,
    adj_p_values_matrix      = adj_p_values_matrix,
    fgsea_results_list       = fgsea_results_list,
    tidy_df                  = tidy_df
  )
}

# ============================================================================
# HEATMAP HELPERS (ATAC STYLE)
# ============================================================================

create_enrichment_plot <- function(enrichment_scores_matrix, p_values_matrix,
                                   adj_p_values_matrix, min_score_global,
                                   max_score_global, title_text,
                                   gwas_gene_counts, deg_counts,
                                   gwas_xlim = NULL, deg_xlim = NULL) {
  enrichment_scores_matrix[is.na(enrichment_scores_matrix)] <- 0
  p_values_matrix[is.na(p_values_matrix)]                   <- 1
  adj_p_values_matrix[is.na(adj_p_values_matrix)]           <- 1
  
  reordered    <- reorder_disease_rows(enrichment_scores_matrix)
  p_reordered  <- reorder_disease_rows(p_values_matrix)
  ap_reordered <- reorder_disease_rows(adj_p_values_matrix)
  
  enrichment_scores_matrix <- reordered$mat
  p_values_matrix          <- p_reordered$mat
  adj_p_values_matrix      <- ap_reordered$mat
  existing_rows            <- reordered$existing_rows
  label_map                <- reordered$label_map
  row_split                <- reordered$row_split
  
  col_fun1 <- colorRamp2(c(min_score_global, 0, max_score_global), c("white", "white", "red"))
  
  gwas_gene_counts <- gwas_gene_counts[names(gwas_gene_counts) %in% names(label_map)]
  gwas_gene_counts <- gwas_gene_counts[existing_rows]
  names(gwas_gene_counts) <- label_map[existing_rows]
  
  deg_counts <- deg_counts[colnames(enrichment_scores_matrix)]
  
  row_anno <- rowAnnotation(
    "GWAS genes" = anno_barplot(
      gwas_gene_counts,
      gp         = gpar(fill = "sienna4"),
      border     = FALSE,
      bar_width  = 0.8,
      axis_param = list(side = "top", labels_rot = 0,
                        at = if (!is.null(gwas_xlim)) seq(0, gwas_xlim[2], length.out = 5) else NULL),
      width = unit(2, "cm"),
      ylim  = gwas_xlim
    )
  )
  
  col_anno <- HeatmapAnnotation(
    "DEG count" = anno_barplot(
      deg_counts,
      gp         = gpar(fill = "sienna4"),
      border     = FALSE,
      bar_width  = 0.8,
      axis_param = list(side = "right",
                        at = if (!is.null(deg_xlim)) seq(0, deg_xlim[2], length.out = 5) else NULL),
      height = unit(2, "cm"),
      ylim   = deg_xlim
    )
  )
  
  Heatmap(
    enrichment_scores_matrix,
    col              = col_fun1,
    top_annotation   = col_anno,
    left_annotation  = row_anno,
    heatmap_legend_param = list(title = "Enrichment\nscore"),
    rect_gp          = gpar(fill = "white", col = "grey90"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      score      <- enrichment_scores_matrix[i, j]
      pval       <- p_values_matrix[i, j]
      adjp       <- adj_p_values_matrix[i, j]
      sq_size    <- ifelse(pval < 0.05, 0.9, 0.7)
      border_col <- ifelse(pval < 0.05, "black", "white")
      grid.rect(x, y, width = width * sq_size, height = height * sq_size,
                gp = gpar(fill = col_fun1(score), col = border_col, lwd = 1.5))
      if (adjp < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 14))
    },
    row_split       = row_split,
    row_gap         = unit(3, "mm"),
    row_title_gp    = gpar(fontsize = 8, fontface = "bold", col = "grey40"),
    row_title_rot   = 0,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title    = title_text,
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  )
}

create_prenatal_postnatal_plot <- function(enrichment_scores_matrix, p_values_matrix,
                                           adj_p_values_matrix, min_score_global,
                                           max_score_global, title_text,
                                           gwas_gene_counts = NULL,
                                           deg_counts = NULL,
                                           gwas_xlim = NULL,
                                           deg_xlim = NULL) {
  enrichment_scores_matrix[is.na(enrichment_scores_matrix)] <- 0
  p_values_matrix[is.na(p_values_matrix)]                   <- 1
  adj_p_values_matrix[is.na(adj_p_values_matrix)]           <- 1
  
  reordered    <- reorder_disease_rows(enrichment_scores_matrix)
  p_reordered  <- reorder_disease_rows(p_values_matrix)
  ap_reordered <- reorder_disease_rows(adj_p_values_matrix)
  
  enrichment_scores_matrix <- reordered$mat
  p_values_matrix          <- p_reordered$mat
  adj_p_values_matrix      <- ap_reordered$mat
  existing_rows            <- reordered$existing_rows
  label_map                <- reordered$label_map
  row_split                <- reordered$row_split
  
  col_fun1 <- colorRamp2(c(min_score_global, 0, max_score_global), c("white", "pink", "red"))
  
  left_anno <- NULL
  if (!is.null(gwas_gene_counts)) {
    gwas_counts_ordered <- gwas_gene_counts[existing_rows]
    names(gwas_counts_ordered) <- label_map[existing_rows]
    if (is.null(gwas_xlim)) gwas_xlim <- c(0, max(gwas_counts_ordered, na.rm = TRUE) * 1.05)
    left_anno <- rowAnnotation(
      `GWAS genes` = anno_barplot(
        gwas_counts_ordered,
        gp         = gpar(fill = "sienna4"),
        border     = FALSE,
        bar_width  = 0.8,
        axis_param = list(side = "top", labels_rot = 0),
        width = unit(2, "cm"),
        ylim  = gwas_xlim
      ),
      show_annotation_name = TRUE,
      annotation_name_gp   = gpar(fontsize = 10),
      annotation_name_rot  = 0
    )
  }
  
  top_anno <- NULL
  if (!is.null(deg_counts)) {
    stages             <- colnames(enrichment_scores_matrix)
    deg_counts_ordered <- deg_counts[stages]
    if (is.null(deg_xlim)) deg_xlim <- c(0, max(deg_counts_ordered, na.rm = TRUE) * 1.05)
    top_anno <- HeatmapAnnotation(
      `DEG count` = anno_barplot(
        deg_counts_ordered,
        gp         = gpar(fill = "sienna4"),
        border     = FALSE,
        bar_width  = 0.8,
        axis_param = list(side = "right", labels_rot = 0),
        height = unit(2, "cm"),
        ylim   = deg_xlim
      ),
      show_annotation_name = TRUE,
      annotation_name_gp   = gpar(fontsize = 10),
      annotation_name_rot  = 0
    )
  }
  
  Heatmap(
    enrichment_scores_matrix,
    col              = col_fun1,
    top_annotation   = top_anno,
    left_annotation  = left_anno,
    heatmap_legend_param = list(title = "Enrichment\nscore"),
    rect_gp          = gpar(fill = "white", col = "grey90"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      score  <- enrichment_scores_matrix[i, j]
      pval   <- p_values_matrix[i, j]
      adjp   <- adj_p_values_matrix[i, j]
      sq_size    <- ifelse(pval < 0.05, 0.9, 0.7)
      border_col <- ifelse(pval < 0.05, "black", "white")
      grid.rect(x, y, width = width * sq_size, height = height * sq_size,
                gp = gpar(fill = col_fun1(score), col = border_col, lwd = 1.5))
      if (adjp < 0.05) grid.text("*", x, y, gp = gpar(fontsize = 14))
    },
    row_split       = row_split,
    row_gap         = unit(3, "mm"),
    row_title_gp    = gpar(fontsize = 8, fontface = "bold", col = "grey40"),
    row_title_rot   = 0,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    row_names_side  = "left",
    row_names_gp    = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title    = title_text,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    width           = unit(2, "cm"),
    column_names_rot = 45,
    show_heatmap_legend = TRUE
  )
}

# ============================================================================
# PIPELINE 1: RNA PRENATAL vs POSTNATAL (ATAC-LIKE INTERLEAVED HEATMAP)
# ============================================================================

run_enrichment_pipeline <- function(prenatal_file, postnatal_file, cell_types,
                                    output_file, pct_threshold = 0.1,
                                    label = "enrichment", interleave = TRUE) {
  output_csv_dir <- file.path("enrichment_csv_results", label)
  
  prenatal_deg  <- prepare_deg_data(read_csv(prenatal_file),  pct_threshold)
  postnatal_deg <- prepare_deg_data(read_csv(postnatal_file), pct_threshold)
  
  prenatal_deg_list  <- select_per_cell(prenatal_deg,  cell_types)
  postnatal_deg_list <- select_per_cell(postnatal_deg, cell_types)
  
  diseases <- c("ASD", "ADHD", "DD", "TS", "GEPI", "CEPI", "SCZ", "MDD", "AN", "OCD",
                "PANIC", "BIP", "PTSD", "PAU", "SUD", "ALS", "FTLD", "MS", "PD", "ALZ",
                "IQ", "EA", "AOW", "Height")
  
  disease_gene <- setNames(lapply(diseases, function(d) toupper(load_disease_genes(d))), diseases)
  
  ASD_rare <- read.csv("rare_variant/ASD_rare.csv")
  ASD_rare <- ASD_rare[(ASD_rare$ASD72 == TRUE | ASD_rare$ASD185 == TRUE) &
                         ASD_rare$FDR_TADA_ASD <= 0.05, ]$gene
  NDD_rare <- read.csv("rare_variant/NDD_rare.csv")
  NDD_rare <- NDD_rare[NDD_rare$significant == TRUE, ]$symbol
  disease_gene[["ASD_rare"]] <- ASD_rare
  disease_gene[["NDD_rare"]] <- NDD_rare
  
  gwas_counts <- sapply(disease_gene, length)
  
  deg_counts_prenatal  <- sapply(prenatal_deg_list,  function(x) length(unique(x$gene)))
  deg_counts_postnatal <- sapply(postnatal_deg_list, function(x) length(unique(x$gene)))
  
  prenatal_ranked <- lapply(cell_types, function(ct) select_significant_degs(prenatal_deg_list[[ct]]))
  names(prenatal_ranked) <- cell_types
  postnatal_ranked <- lapply(cell_types, function(ct) select_significant_degs(postnatal_deg_list[[ct]]))
  names(postnatal_ranked) <- cell_types
  
  prenatal_results  <- run_fgsea(cell_types, disease_gene, prenatal_ranked,  'prenatal',
                                 label = label, output_csv_dir = output_csv_dir)
  postnatal_results <- run_fgsea(cell_types, disease_gene, postnatal_ranked, 'postnatal',
                                 label = label, output_csv_dir = output_csv_dir)
  
  # Excel output
  wb <- createWorkbook()
  all_tidy <- bind_rows(prenatal_results$tidy_df, postnatal_results$tidy_df)
  append_sheet(wb, all_tidy,                  sheet_name = "All_results")
  append_sheet(wb, prenatal_results$tidy_df,   sheet_name = "Prenatal")
  append_sheet(wb, postnatal_results$tidy_df,  sheet_name = "Postnatal")
  excel_file <- sub("\\.pdf$", "_results.xlsx", output_file)
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  all_scores <- c(as.vector(prenatal_results$enrichment_scores_matrix),
                  as.vector(postnatal_results$enrichment_scores_matrix))
  min_score_global <- min(all_scores, na.rm = TRUE)
  max_score_global <- max(all_scores, na.rm = TRUE)
  
  legend_pval <- Legend(
    title     = "p-value",
    labels    = c("< 0.05", "≥ 0.05"),
    type      = "points",
    pch       = 15,
    size      = unit(c(6, 3), "mm"),
    legend_gp = gpar(fill = c("grey", "grey"), col = 'grey')
  )
  
  legend_star <- Legend(
    title     = "Adjusted p-value",
    labels    = c("adj(p) < 0.05"),
    type      = "points",
    legend_gp = gpar(fontsize = 10),
    labels_gp = gpar(col = "black")
  )
  
  if (interleave) {
    nes_int  <- interleave_pre_post_matrices(prenatal_results$enrichment_scores_matrix,
                                             postnatal_results$enrichment_scores_matrix,
                                             cell_types)
    pval_int <- interleave_pre_post_matrices(prenatal_results$p_values_matrix,
                                             postnatal_results$p_values_matrix,
                                             cell_types)
    padj_int <- interleave_pre_post_matrices(prenatal_results$adj_p_values_matrix,
                                             postnatal_results$adj_p_values_matrix,
                                             cell_types)
    deg_int  <- interleave_pre_post_vector(deg_counts_prenatal, deg_counts_postnatal, cell_types)
    
    ht <- create_enrichment_plot(
      nes_int, pval_int, padj_int,
      min_score_global, max_score_global,
      title_text = label,
      gwas_gene_counts = gwas_counts,
      deg_counts       = deg_int
    )
    
    pdf_w <- max(12, length(cell_types) * 1.5 + 5)
    pdf(output_file, width = pdf_w, height = 14)
    draw(ht,
         annotation_legend_list = list(legend_pval, legend_star),
         heatmap_legend_side    = "right",
         annotation_legend_side = "bottom")
    dev.off()
  } else {
    ht1 <- create_enrichment_plot(
      prenatal_results$enrichment_scores_matrix,
      prenatal_results$p_values_matrix,
      prenatal_results$adj_p_values_matrix,
      min_score_global, max_score_global, "Prenatal",
      gwas_gene_counts = gwas_counts,
      deg_counts       = deg_counts_prenatal
    )
    ht2 <- create_enrichment_plot(
      postnatal_results$enrichment_scores_matrix,
      postnatal_results$p_values_matrix,
      postnatal_results$adj_p_values_matrix,
      min_score_global, max_score_global, "Postnatal",
      gwas_gene_counts = gwas_counts,
      deg_counts       = deg_counts_postnatal
    )
    
    pdf_w <- max(13, length(cell_types) * 1.1 + 8)
    pdf(output_file, width = pdf_w, height = 14)
    draw(ht1 + ht2,
         annotation_legend_list = list(legend_pval, legend_star),
         heatmap_legend_side    = "right",
         annotation_legend_side = "bottom")
    dev.off()
  }
  
  cat("\n========================================\n")
  cat("Analysis complete!\n")
  cat("PDF saved as:        ", output_file, "\n")
  cat("Excel saved as:      ", excel_file, "\n")
  cat("CSV results saved in:", output_csv_dir, "\n")
  cat("========================================\n\n")
  
  list(prenatal = prenatal_results, postnatal = postnatal_results)
}

# ============================================================================
# PIPELINE 2: PRENATAL vs POSTNATAL BY CELL TYPE (row_split + row_gap)
# ============================================================================

run_prenatal_postnatal_pipeline <- function(marker_directory, output_file,
                                            selected_cell_types = NULL,
                                            pct_threshold = 0.1,
                                            label = "prenatal_postnatal") {
  cat("\n========================================\n")
  cat("Prenatal vs Postnatal Enrichment Analysis\n")
  cat("========================================\n\n")
  
  output_csv_dir <- file.path("enrichment_csv_results", label)
  
  cat("Loading marker genes from", marker_directory, "...\n")
  marker_data <- load_marker_genes(marker_directory)
  colnames(marker_data) <- c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2",
                             "p_val_adj", "Stage", "term")
  
  split_term <- do.call(rbind, strsplit(marker_data$term, "_"))
  marker_data$cell_type <- split_term[, 1]
  marker_data <- marker_data %>% filter(Stage %in% c("Prenatal", "Postnatal"))
  
  marker_data$cell_type <- recode(marker_data$cell_type,
                                  "Micro" = "Microglia",
                                  "Sst"   = "SST",
                                  "Vip"   = "VIP",
                                  "Pvalb" = "PVALB")
  
  marker_data$Stage <- factor(marker_data$Stage, levels = c("Postnatal", "Prenatal"))
  
  available_cell_types <- unique(marker_data$cell_type)
  cat("Available cell types:", paste(available_cell_types, collapse = ", "), "\n\n")
  
  if (is.null(selected_cell_types)) {
    if (interactive()) {
      cat("Available: ", paste(available_cell_types, collapse = ", "), "\n")
      user_input <- readline(prompt = "> ")
      if (nchar(trimws(user_input)) == 0) {
        selected_cell_types <- available_cell_types
      } else {
        selected_cell_types <- trimws(unlist(strsplit(user_input, ",")))
      }
    } else {
      selected_cell_types <- available_cell_types
    }
  }
  
  invalid_types <- setdiff(selected_cell_types, available_cell_types)
  if (length(invalid_types) > 0) {
    stop("Invalid cell types: ", paste(invalid_types, collapse = ", "))
  }
  
  cell_types <- selected_cell_types
  cat("Analyzing cell types:", paste(cell_types, collapse = ", "), "\n")
  cat("Using pct_threshold:", pct_threshold, "\n\n")
  
  diseases <- c("ASD", "ADHD", "DD", "TS", "GEPI", "CEPI", "SCZ", "MDD", "AN", "OCD",
                "PANIC", "BIP", "PTSD", "PAU", "SUD", "ALS", "FTLD", "MS", "PD", "ALZ",
                "IQ", "EA", "AOW", "Height")
  
  disease_gene <- setNames(lapply(diseases, function(disease) {
    cat("  Loading", disease, "...\n")
    toupper(load_disease_genes(disease))
  }), diseases)
  
  cat("  Loading rare variant genes...\n")
  ASD_rare <- read.csv("rare_variant/ASD_rare.csv", stringsAsFactors = FALSE)
  ASD_rare <- ASD_rare[(ASD_rare$ASD72 == TRUE | ASD_rare$ASD185 == TRUE) &
                         ASD_rare$FDR_TADA_ASD <= 0.05, ]$gene
  NDD_rare <- read.csv("rare_variant/NDD_rare.csv", stringsAsFactors = FALSE)
  NDD_rare <- NDD_rare[NDD_rare$significant == TRUE, ]$symbol
  disease_gene[["ASD_rare"]] <- ASD_rare
  disease_gene[["NDD_rare"]] <- NDD_rare
  
  gwas_counts <- sapply(disease_gene, length)
  
  cat("\nFiltering genes by pct_threshold (", pct_threshold, ")...\n", sep = "")
  marker_data_filtered <- marker_data %>% filter(pct.1 >= pct_threshold | pct.2 >= pct_threshold)
  
  results_list <- list()
  for (cell_type in cell_types) {
    cat("\n=== Analyzing", cell_type, "===\n")
    results_list[[cell_type]] <- run_fgsea_prenatal_postnatal(
      marker_data    = marker_data_filtered,
      disease_gene   = disease_gene,
      cell_type      = cell_type,
      label          = label,
      output_csv_dir = output_csv_dir
    )
  }
  
  # Excel output
  wb <- createWorkbook()
  all_tidy <- bind_rows(lapply(results_list, function(x) x$tidy_df))
  append_sheet(wb, all_tidy, sheet_name = "All_results")
  for (ct in names(results_list)) {
    ct_tidy <- results_list[[ct]]$tidy_df
    if (nrow(ct_tidy) > 0) append_sheet(wb, ct_tidy, sheet_name = substr(ct, 1, 31))
  }
  excel_file <- sub("\\.pdf$", "_results.xlsx", output_file)
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("\nExcel saved as:", excel_file, "\n")
  
  all_scores <- unlist(lapply(results_list, function(x) as.vector(x$enrichment_scores_matrix)))
  all_scores <- all_scores[!is.na(all_scores)]
  min_score_global <- min(all_scores)
  max_score_global <- max(all_scores)
  
  cat("\nCalculating DEG counts per cell type and stage...\n")
  deg_counts_per_cell <- list()
  for (cell_type in cell_types) {
    cell_markers <- marker_data_filtered %>% filter(cell_type == !!cell_type)
    deg_counts_per_cell[[cell_type]] <- c(
      "Postnatal" = length(unique(cell_markers[cell_markers$Stage == "Postnatal", ]$gene)),
      "Prenatal"  = length(unique(cell_markers[cell_markers$Stage == "Prenatal",  ]$gene))
    )
  }
  
  all_deg_counts <- unlist(deg_counts_per_cell)
  deg_xlim  <- c(0, max(all_deg_counts, na.rm = TRUE) * 1.05)
  gwas_xlim <- c(0, max(gwas_counts,    na.rm = TRUE) * 1.05)
  
  cat("Creating heatmaps...\n")
  heatmap_list <- list()
  for (cell_type in cell_types) {
    cat("  Creating heatmap for", cell_type, "...\n")
    heatmap_list[[cell_type]] <- create_prenatal_postnatal_plot(
      enrichment_scores_matrix = results_list[[cell_type]]$enrichment_scores_matrix,
      p_values_matrix          = results_list[[cell_type]]$p_values_matrix,
      adj_p_values_matrix      = results_list[[cell_type]]$adj_p_values_matrix,
      min_score_global         = min_score_global,
      max_score_global         = max_score_global,
      title_text               = cell_type,
      gwas_gene_counts         = gwas_counts,
      deg_counts               = deg_counts_per_cell[[cell_type]],
      gwas_xlim                = gwas_xlim,
      deg_xlim                 = deg_xlim
    )
  }
  
  legend_pval <- Legend(
    title     = "p-value",
    labels    = c("< 0.05", "≥ 0.05"),
    type      = "points",
    pch       = 15,
    size      = unit(c(6, 3), "mm"),
    legend_gp = gpar(fill = c("grey", "grey"), col = 'grey')
  )
  legend_star <- Legend(
    title     = "Adjusted p-value",
    labels    = c("adj(p) < 0.05"),
    type      = "points",
    legend_gp = gpar(fontsize = 10),
    labels_gp = gpar(col = "black")
  )
  
  cat("\nSaving plot to", output_file, "...\n")
  pdf(output_file, width = 20, height = 6)
  combined_ht <- heatmap_list[[1]]
  if (length(heatmap_list) > 1) {
    for (i in 2:length(heatmap_list)) combined_ht <- combined_ht + heatmap_list[[i]]
  }
  draw(combined_ht,
       annotation_legend_list = list(legend_pval, legend_star),
       heatmap_legend_side    = "right",
       annotation_legend_side = "bottom")
  dev.off()
  
  cat("\n========================================\n")
  cat("Analysis complete!\n")
  cat("PDF saved as:", output_file, "\n")
  cat("Excel saved as:", excel_file, "\n")
  cat("CSV results in:", output_csv_dir, "\n")
  cat("========================================\n\n")
  
  results_list
}

# ============================================================================
# RUN
# ============================================================================

run_enrichment_pipeline(
  prenatal_file  = "RNAmarkers1692025/prenatal_markers_neuronvsnn_celltypes_all_sep24.csv",
  postnatal_file = "RNAmarkers1692025/postnatal_markers_neuronvsnn_celltypes_all_sep24.csv",
  cell_types     = c("Neuron", "Micro", "OPC", "Astrocyte"),
  output_file    = "RNA_enrichment_heatmaps_neuron_nonneuron_correct_counts.pdf",
  pct_threshold  = 0.2,
  label          = "neuronvsnn_sep24",
  interleave     = TRUE
)

run_enrichment_pipeline(
  prenatal_file  = "RNAmarkers1692025/prenatal_markers_minor_celltypes_all.csv",
  postnatal_file = "RNAmarkers1692025/postnatal_markers_minor_celltypes_all.csv",
  cell_types     = c(
    "L6 CT", "L4 IT", "L5/6 NP", "L5 IT", "L6b", "L2/3 IT",
    "L6 IT Car3", "L6 IT", "L5 ET", "Sncg", "Lamp5", "Vip",
    "Chandelier", "Pvalb", "Sst", "Lamp5 Lhx6", "Pax6", "Sst Chodl",
    "Astrocyte", "Micro", "OPC", "PC", "Oligo", "VLMC", "Immune", "Endo", "SMC"
  ),
  output_file   = "RNA_enrichment_heatmaps_cell_subtype_correct_counts.pdf",
  pct_threshold = 0.2,
  label         = "minor_celltypes",
  interleave    = FALSE
)

run_enrichment_pipeline(
  prenatal_file  = "RNAmarkers1692025/markers_dev_cell_prenatal.csv",
  postnatal_file = "RNAmarkers1692025/markers_dev_cell_postnatal.csv",
  cell_types     = c("ExN", "InN", "Micro", "OPC", "Astrocyte"),
  output_file    = "RNA_enrichment_heatmaps_dev_cell_correct_counts.pdf",
  pct_threshold  = 0.2,
  label          = "dev_cell",
  interleave     = TRUE
)

prenatal_postnatal_results <- run_prenatal_postnatal_pipeline(
  marker_directory    = "markers_9152025/",
  output_file         = "RNA_enrichment_heatmaps_prenatal_vs_postnatal.pdf",
  selected_cell_types = c("ExN", "InN", "Glia", "NNC", "SST", "VIP", "PVALB",
                          "Microglia", "OPC", "PC", "VLMC"),
  pct_threshold = 0.2,
  label         = "markers_9152025_prenatal_postnatal"
)
