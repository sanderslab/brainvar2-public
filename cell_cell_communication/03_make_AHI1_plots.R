#!/usr/bin/env Rscript

# ============================================================
# 03_make_AHI1_ALL_plots_FINAL_v3.R
#
# ONE downstream plotting script for Steps 1 + 2.
#
# It DOES NOT rerun CellChat inference.
# It reads:
#   1) per-sample CellChat RDS objects from Step 1
#   2) AHI1_plot_inputs.rds + differential tables from Step 2
#
# Main outputs:
#   01_differential/     volcanoes, labeled volcanoes, heatmaps
#   02_lr_boxplots/      top LR healthy-vs-risk plots
#   03_pathway_plots/    pathway volcanoes, heatmaps, boxplots
#   04_networks/         per-sample count/weight networks
#   05_cellchat_paths/   standard CellChat pathway visualizations
#   06_group_networks/   healthy/risk averaged network comparisons
#   07_targeted/         glutamate/cholinergic + ExN/InN -> OPC plots
#   08_tables/            copied/derived result tables
#   09_manifest/          list of generated files and settings
#
# The script is deliberately defensive: unavailable CellChat plots are
# skipped with a warning rather than stopping the complete plotting run.
# ============================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(grid)
})

options(stringsAsFactors = FALSE)

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------
BASE_OUT_DIR <- "cellchat_results_01"
MODEL_NAME <- "truncatedMean10"
GROUP_DEFINITION <- "subtype4"
RESTRICT_STAGE <- "prenatal"

MIN_PROBABILITY <- 0.001
MAX_TOP_LR <- 20
MAX_TOP_PATHWAYS <- 20
MAX_PATHWAYS_FOR_PER_SAMPLE <- 12
MAX_SAMPLES_FOR_FULL_PER_SAMPLE_PLOTS <- Inf

MAKE_PER_SAMPLE_PLOTS <- TRUE
MAKE_STANDARD_CELLCHAT_PLOTS <- TRUE
MAKE_GROUP_AVERAGE_NETWORKS <- TRUE
MAKE_TARGETED_PLOTS <- TRUE

# Pathways to prioritize. Additional detected pathways are also plotted.
TARGET_PATHWAYS <- c("Glutamate", "Ach")

# Targeted cell types for the H1 network view.
H1_SENDERS <- c("ExN", "InN")
H1_RECEIVER <- "OPC"

results_dir <- file.path(
  BASE_OUT_DIR,
  "homo_differential",
  MODEL_NAME,
  GROUP_DEFINITION,
  if (is.null(RESTRICT_STAGE)) "all_stages" else RESTRICT_STAGE
)

cellchat_dir <- file.path(
  BASE_OUT_DIR,
  "cellchat_objects",
  MODEL_NAME,
  GROUP_DEFINITION
)

plot_root <- file.path(results_dir, "ALL_PLOTS")
plot_dirs <- list(
  differential = file.path(plot_root, "01_differential"),
  lr_boxplots = file.path(plot_root, "02_lr_boxplots"),
  pathway_plots = file.path(plot_root, "03_pathway_plots"),
  networks = file.path(plot_root, "04_networks"),
  cellchat_paths = file.path(plot_root, "05_cellchat_paths"),
  group_networks = file.path(plot_root, "06_group_networks"),
  targeted = file.path(plot_root, "07_targeted"),
  tables = file.path(plot_root, "08_tables"),
  manifest = file.path(plot_root, "09_manifest")
)

if (!is.list(plot_dirs)) stop("INTERNAL ERROR: plot_dirs must be a list. Please use this FINAL_v2 script.")
if (is.null(names(plot_dirs))) stop("INTERNAL ERROR: plot_dirs has no names.")
for (d in unname(plot_dirs)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

msg <- function(...) cat(paste0(..., "\n"))
msg("RUNNING SCRIPT: 03_make_AHI1_ALL_plots_FINAL_v4.R")

safe_write_csv <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  tryCatch({
    readr::write_csv(as.data.frame(x), file)
    TRUE
  }, error = function(e) {
    warning("Could not write ", file, ": ", conditionMessage(e))
    FALSE
  })
}

safe_ggsave <- function(file, plot, width, height, ...) {
  tryCatch({
    ggplot2::ggsave(file, plot, width = width, height = height, ...)
    TRUE
  }, error = function(e) {
    warning("Plot failed: ", basename(file), " -> ", conditionMessage(e))
    FALSE
  })
}

safe_pdf_call <- function(file, expr, width = 7, height = 7) {
  opened <- FALSE
  ok <- FALSE
  tryCatch({
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    opened <- TRUE
    force(expr)
    ok <- TRUE
  }, error = function(e) {
    warning("Plot failed: ", basename(file), " -> ", conditionMessage(e))
  }, finally = {
    if (opened && grDevices::dev.cur() > 1) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  })
  ok
}

safe_cellchat_plot <- function(file, FUN, width = 7, height = 7) {
  opened <- FALSE
  ok <- FALSE
  tryCatch({
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    opened <- TRUE
    FUN()
    ok <- TRUE
  }, error = function(e) {
    warning("CellChat plot failed: ", basename(file), " -> ", conditionMessage(e))
  }, finally = {
    if (opened && grDevices::dev.cur() > 1) {
      try(grDevices::dev.off(), silent = TRUE)
    }
  })
  ok
}

sanitize <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  ifelse(nzchar(x), x, "unknown")
}

add_sig_class <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df %>%
    mutate(
      P.Value = as.numeric(P.Value),
      adj.P.Val = as.numeric(adj.P.Val),
      logFC = as.numeric(logFC),
      neglog10P = -log10(pmax(P.Value, .Machine$double.xmin)),
      significance = case_when(
        is.finite(adj.P.Val) & adj.P.Val < 0.05 ~ "FDR < 0.05",
        is.finite(P.Value) & P.Value < 0.05 ~ "P < 0.05",
        TRUE ~ "ns"
      )
    )
}

prepare_heatmap <- function(mat, max_rows = 20) {
  if (is.null(mat) || length(mat) == 0 || nrow(mat) == 0 || ncol(mat) == 0) return(NULL)
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  keep <- apply(mat, 1, function(x) {
    x <- x[is.finite(x)]
    length(x) >= 2 && is.finite(sd(x)) && sd(x) > 0
  })
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)
  if (nrow(mat) > max_rows) mat <- mat[seq_len(max_rows), , drop = FALSE]
  z <- t(scale(t(mat)))
  z[!is.finite(z)] <- 0
  z
}

annotation_for <- function(coldata, samples) {
  x <- as.data.frame(coldata)
  if (!"Sample" %in% colnames(x)) return(NULL)
  x$Sample <- trimws(as.character(x$Sample))
  x <- x %>% distinct(Sample, .keep_all = TRUE) %>% filter(Sample %in% samples)
  x <- x[match(samples, x$Sample), , drop = FALSE]
  rownames(x) <- x$Sample
  keep <- intersect(c("AHI1_status", "stage"), colnames(x))
  if (length(keep) == 0) return(NULL)
  x[, keep, drop = FALSE]
}

p_label <- function(p) {
  if (is.na(p) || !is.finite(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", formatC(p, format = "f", digits = 3))
}

wilcox_safe <- function(y, group) {
  ok <- is.finite(y) & !is.na(group)
  y <- y[ok]
  group <- factor(group[ok], levels = c("healthy", "risk"))
  if (sum(group == "healthy") < 2 || sum(group == "risk") < 2) return(NA_real_)
  if (all(y[group == "healthy"] == y[group == "risk"])) return(1)
  tryCatch(wilcox.test(y ~ group, exact = FALSE)$p.value, error = function(e) NA_real_)
}

make_boxplot <- function(df, value_col, title, ylab, p_value = NA_real_) {
  d <- df %>%
    transmute(
      Sample = Sample,
      AHI1_status = factor(AHI1_status, levels = c("healthy", "risk")),
      value = as.numeric(.data[[value_col]])
    ) %>%
    filter(!is.na(AHI1_status), is.finite(value))
  if (nrow(d) == 0) return(NULL)
  ymax <- max(d$value, na.rm = TRUE)
  ymin <- min(d$value, na.rm = TRUE)
  span <- ymax - ymin
  if (!is.finite(span) || span == 0) span <- max(abs(ymax), 1) * 0.1
  ggplot(d, aes(AHI1_status, value, fill = AHI1_status)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.75) +
    annotate("text", x = 1.5, y = ymax + 0.12 * span, label = p_label(p_value), size = 4) +
    scale_fill_manual(values = c(healthy = "#1b9e77", risk = "#d95f02")) +
    coord_cartesian(ylim = c(ymin, ymax + 0.25 * span), clip = "off") +
    labs(title = title, x = NULL, y = ylab) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
}

# ------------------------------------------------------------
# LOAD STEP-2 PLOT INPUTS
# ------------------------------------------------------------

bundle_file <- file.path(results_dir, "AHI1_plot_inputs.rds")
if (!file.exists(bundle_file)) {
  stop("Missing: ", bundle_file, ". Run Step 2 first.")
}

bundle <- readRDS(bundle_file)
cci_mat <- bundle$cci_mat
path_mat <- bundle$path_mat
coldata <- as.data.frame(bundle$coldata)
lr_res <- bundle$lr_results
path_res <- bundle$pathway_results
path_global_res <- bundle$pathway_global_results

if ("Sample" %in% colnames(coldata)) coldata$Sample <- trimws(as.character(coldata$Sample))
if (!is.null(cci_mat)) cci_mat <- as.matrix(cci_mat)
if (!is.null(path_mat)) path_mat <- as.matrix(path_mat)

msg("Loaded Step-2 plot inputs.")
msg("Results directory: ", results_dir)
msg("Plot directory:    ", plot_root)

# ------------------------------------------------------------
# LOAD PER-SAMPLE CELLCHAT OBJECTS FROM STEP 1
# ------------------------------------------------------------

rds_files <- list.files(cellchat_dir, pattern = "\\.rds$", full.names = TRUE)
if (length(rds_files) == 0) stop("No Step-1 CellChat RDS files found in: ", cellchat_dir)

obj_list <- list()
for (f in rds_files) {
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj) || !inherits(obj, "CellChat")) next
  sid <- NA_character_
  if (!is.null(obj@meta$samples)) {
    s <- unique(as.character(obj@meta$samples))
    if (length(s) == 1) sid <- s
  }
  if (is.na(sid)) sid <- sub("^CCI_CellChat_.*_", "", tools::file_path_sans_ext(basename(f)))
  obj_list[[sid]] <- obj
}

if (length(obj_list) == 0) stop("No readable CellChat objects found.")

# Keep samples that are represented in Step 2, when available.
if ("Sample" %in% colnames(coldata)) {
  keep_names <- intersect(names(obj_list), unique(coldata$Sample))
  if (length(keep_names) > 0) obj_list <- obj_list[keep_names]
}

msg("Loaded ", length(obj_list), " CellChat objects from Step 1.")

# ------------------------------------------------------------
# 00. MANIFEST / SUMMARY
# ------------------------------------------------------------

sample_summary <- coldata %>%
  count(AHI1_status, name = "n_samples")

safe_write_csv(sample_summary, file.path(plot_dirs$tables, "00_sample_counts.csv"))

object_summary <- bind_rows(lapply(names(obj_list), function(sid) {
  x <- obj_list[[sid]]
  g <- tryCatch(as.character(x@idents), error = function(e) character(0))
  data.frame(
    Sample = sid,
    n_cells = ncol(x@data),
    groups = paste(sort(unique(g)), collapse = ";"),
    n_groups = length(unique(g)),
    n_LR = if (!is.null(x@net$prob)) dim(x@net$prob)[3] else 0,
    n_pathways = if (!is.null(x@netP$prob)) dim(x@netP$prob)[3] else 0
  )
}))
safe_write_csv(object_summary, file.path(plot_dirs$tables, "00_CellChat_object_summary.csv"))

# ------------------------------------------------------------
# 01. DIFFERENTIAL LR PLOTS
# ------------------------------------------------------------

if (!is.null(lr_res) && nrow(lr_res) > 0) {
  lr_tbl <- add_sig_class(lr_res)
  safe_write_csv(lr_tbl, file.path(plot_dirs$tables, "01_LR_differential_results.csv"))

  p_lr <- ggplot(lr_tbl, aes(logFC, neglog10P)) +
    geom_point(aes(color = significance), alpha = 0.75, size = 1.5, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c(ns = "grey70", `P < 0.05` = "#1f78b4", `FDR < 0.05` = "#e31a1c")) +
    theme_bw(base_size = 11) +
    labs(title = "AHI1 differential CellChat LR settings", x = "Risk - healthy coefficient", y = "-log10(P)", color = NULL)

  safe_ggsave(file.path(plot_dirs$differential, "01_volcano_LR.pdf"), p_lr, 8, 6)

  top_lr <- lr_tbl %>% arrange(P.Value) %>% slice_head(n = min(MAX_TOP_LR, nrow(lr_tbl)))
  if (nrow(top_lr) > 0) {
    p_labeled <- p_lr + geom_text_repel(data = top_lr, aes(label = settings), size = 2.6, max.overlaps = Inf, na.rm = TRUE)
    safe_ggsave(file.path(plot_dirs$differential, "02_volcano_LR_labeled.pdf"), p_labeled, 11, 8)

    if (!is.null(cci_mat) && nrow(cci_mat) > 0) {
      ids <- intersect(top_lr$settings, rownames(cci_mat))
      hm <- prepare_heatmap(cci_mat[ids, , drop = FALSE], MAX_TOP_LR)
      if (!is.null(hm)) {
        ann <- annotation_for(coldata, colnames(hm))
        common <- intersect(colnames(hm), rownames(ann))
        if (length(common) > 1) {
          hm <- hm[, common, drop = FALSE]
          ann <- ann[common, , drop = FALSE]
          if ("AHI1_status" %in% colnames(ann)) {
            ord <- order(factor(ann$AHI1_status, levels = c("healthy", "risk")))
            hm <- hm[, ord, drop = FALSE]
            ann <- ann[ord, , drop = FALSE]
          }
          safe_pdf_call(file.path(plot_dirs$differential, "03_heatmap_top20_LR.pdf"), {
            pheatmap(hm, annotation_col = ann, show_colnames = FALSE, cluster_cols = FALSE,
                     clustering_distance_rows = "correlation", main = "Top differential CellChat LR settings")
          }, 12, 9)
        }
      }
    }
  }

  # Top LR individual healthy-vs-risk boxplots.
  if (!is.null(cci_mat) && nrow(cci_mat) > 0 && "Sample" %in% colnames(coldata)) {
    status <- coldata %>% select(Sample, AHI1_status) %>% distinct()
    box_master <- list()
    for (id in top_lr$settings) {
      if (!(id %in% rownames(cci_mat))) next
      d <- data.frame(Sample = colnames(cci_mat), value = as.numeric(cci_mat[id, ])) %>% left_join(status, by = "Sample")
      pval <- wilcox_safe(d$value, d$AHI1_status)
      p <- make_boxplot(d, "value", paste0(id, " - LR communication"), "Communication probability", pval)
      if (!is.null(p)) safe_ggsave(file.path(plot_dirs$lr_boxplots, paste0("LR_boxplot_", sanitize(id), ".pdf")), p, 6, 5.5)
      box_master[[length(box_master) + 1L]] <- data.frame(settings = id, boxplot_P = pval)
    }
    if (length(box_master)) safe_write_csv(bind_rows(box_master), file.path(plot_dirs$tables, "02_top_LR_boxplot_P_values.csv"))
  }
}

# ------------------------------------------------------------
# 02. TARGETED LR TABLE
# ------------------------------------------------------------

if (!is.null(lr_res) && nrow(lr_res) > 0) {
  target_lr <- lr_res %>%
    filter(grepl("Glutamate|Glu-|GRIA|GRIK|GRM|SLC17|SLC1A|GLS|GLS2|Acetylcholine|Ach-|CHAT|SLC5A7|SLC18A3|CHRN|CHRNA|CHRM",
                 settings, ignore.case = TRUE))
  safe_write_csv(target_lr, file.path(plot_dirs$targeted, "01_glutamate_cholinergic_LR.csv"))
}

# ------------------------------------------------------------
# 03. PATHWAY DIFFERENTIAL PLOTS
# ------------------------------------------------------------

if (!is.null(path_res) && nrow(path_res) > 0) {
  path_tbl <- add_sig_class(path_res)
  safe_write_csv(path_tbl, file.path(plot_dirs$tables, "03_pathway_sender_receiver_results.csv"))

  p_path <- ggplot(path_tbl, aes(logFC, neglog10P)) +
    geom_point(aes(color = significance), alpha = 0.75, size = 1.6, na.rm = TRUE) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_color_manual(values = c(ns = "grey70", `P < 0.05` = "#1f78b4", `FDR < 0.05` = "#e31a1c")) +
    theme_bw(base_size = 11) +
    labs(title = "AHI1 differential pathway sender-receiver settings", x = "Risk - healthy coefficient", y = "-log10(P)", color = NULL)

  safe_ggsave(file.path(plot_dirs$pathway_plots, "01_volcano_pathway_sender_receiver.pdf"), p_path, 8, 6)

  top_path_sr <- path_tbl %>% arrange(P.Value) %>% slice_head(n = min(MAX_TOP_PATHWAYS, nrow(path_tbl)))
  if (nrow(top_path_sr) > 0) {
    p_path_lab <- p_path + geom_text_repel(data = top_path_sr, aes(label = settings), size = 2.6, max.overlaps = Inf, na.rm = TRUE)
    safe_ggsave(file.path(plot_dirs$pathway_plots, "02_volcano_pathway_sender_receiver_labeled.pdf"), p_path_lab, 11, 8)

    if (!is.null(path_mat) && nrow(path_mat) > 0) {
      ids <- intersect(top_path_sr$settings, rownames(path_mat))
      hm <- prepare_heatmap(path_mat[ids, , drop = FALSE], MAX_TOP_PATHWAYS)
      if (!is.null(hm)) {
        ann <- annotation_for(coldata, colnames(hm))
        common <- intersect(colnames(hm), rownames(ann))
        if (length(common) > 1) {
          hm <- hm[, common, drop = FALSE]
          ann <- ann[common, , drop = FALSE]
          if ("AHI1_status" %in% colnames(ann)) {
            ord <- order(factor(ann$AHI1_status, levels = c("healthy", "risk")))
            hm <- hm[, ord, drop = FALSE]
            ann <- ann[ord, , drop = FALSE]
          }
          safe_pdf_call(file.path(plot_dirs$pathway_plots, "03_heatmap_top_pathway_sender_receiver.pdf"), {
            pheatmap(hm, annotation_col = ann, show_colnames = FALSE, cluster_cols = FALSE,
                     clustering_distance_rows = "correlation", main = "Top differential pathway sender-receiver settings")
          }, 12, 10)
        }
      }
    }
  }

  target_path <- path_tbl %>% filter(grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline", settings, ignore.case = TRUE))
  safe_write_csv(target_path, file.path(plot_dirs$targeted, "02_glutamate_cholinergic_pathway_sender_receiver.csv"))
}

# ------------------------------------------------------------
# 04. GLOBAL PATHWAY RESULTS + BOXPLOTS
# ------------------------------------------------------------

if (!is.null(path_global_res) && nrow(path_global_res) > 0) {
  global_tbl <- add_sig_class(path_global_res)
  safe_write_csv(global_tbl, file.path(plot_dirs$tables, "04_global_pathway_results.csv"))

  top_global <- global_tbl %>% arrange(P.Value) %>% slice_head(n = min(MAX_TOP_PATHWAYS, nrow(global_tbl)))
  if (nrow(top_global) > 0) {
    p_global <- ggplot(top_global, aes(x = reorder(settings, logFC), y = logFC)) +
      geom_hline(yintercept = 0, linetype = "dashed") + geom_point(size = 2) + coord_flip() +
      theme_bw(base_size = 11) + labs(title = "Top global pathway differences", x = NULL, y = "Risk - healthy coefficient")
    safe_ggsave(file.path(plot_dirs$pathway_plots, "04_top_global_pathways.pdf"), p_global, 9, 8)
  }

  sig_global <- global_tbl %>% filter(is.finite(P.Value), P.Value < 0.05) %>% arrange(P.Value) %>% slice_head(n = MAX_TOP_PATHWAYS)
  safe_write_csv(sig_global, file.path(plot_dirs$pathway_plots, "05_P_significant_global_pathways.csv"))

  target_global <- global_tbl %>% filter(grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline", settings, ignore.case = TRUE))
  safe_write_csv(target_global, file.path(plot_dirs$targeted, "03_glutamate_cholinergic_pathway_global.csv"))

  # Global pathway weight and LR count boxplots, based on the exact raw Step-2 tables.
  raw_path_file <- file.path(results_dir, "pathway_all_samples_raw.csv")
  raw_lr_file <- file.path(results_dir, "LR_all_samples_raw.csv")

  if (file.exists(raw_path_file) && file.exists(raw_lr_file) && nrow(sig_global) > 0) {
    raw_path <- read.csv(raw_path_file, check.names = FALSE)
    raw_lr <- read.csv(raw_lr_file, check.names = FALSE)
    raw_path$value <- as.numeric(raw_path$value)
    raw_lr$value <- as.numeric(raw_lr$value)
    raw_path$Sample <- trimws(as.character(raw_path$Sample))
    raw_lr$Sample <- trimws(as.character(raw_lr$Sample))

    # Map LR -> pathway using the current CellChat database.
    db_int <- tryCatch(CellChatDB.human$interaction %>% select(interaction_name, pathway_name) %>% distinct(), error = function(e) NULL)
    if (!is.null(db_int)) raw_lr <- raw_lr %>% left_join(db_int, by = c("LR" = "interaction_name"))

    status <- coldata %>% select(Sample, AHI1_status) %>% distinct()
    weight_df <- raw_path %>%
      group_by(Pathway, Sample) %>%
      summarise(weight = sum(value, na.rm = TRUE), .groups = "drop")
    count_df <- raw_lr %>%
      filter(!is.na(pathway_name), is.finite(value), value > 0) %>%
      group_by(pathway_name, Sample) %>%
      summarise(n_LR = n_distinct(LR), .groups = "drop")

    master_box <- list()
    for (pwy in sig_global$settings) {
      w <- weight_df %>% filter(Pathway == pwy) %>% right_join(status, by = "Sample") %>% mutate(value = replace_na(weight, 0))
      c <- count_df %>% filter(pathway_name == pwy) %>% right_join(status, by = "Sample") %>% mutate(value = replace_na(n_LR, 0))

      wp <- wilcox_safe(w$value, w$AHI1_status)
      cp <- wilcox_safe(c$value, c$AHI1_status)

      pw <- make_boxplot(w, "value", paste0(pwy, " - pathway communication weight"), "Summed pathway communication probability", wp)
      pc <- make_boxplot(c, "value", paste0(pwy, " - number of LR interactions"), "Number of detected LR interactions", cp)

      s <- sanitize(pwy)
      if (!is.null(pw)) safe_ggsave(file.path(plot_dirs$pathway_plots, paste0("06_boxplot_weight_", s, ".pdf")), pw, 6, 5.5)
      if (!is.null(pc)) safe_ggsave(file.path(plot_dirs$pathway_plots, paste0("07_boxplot_LRcount_", s, ".pdf")), pc, 6, 5.5)

      master_box[[length(master_box) + 1L]] <- data.frame(
        Pathway = pwy,
        global_P = sig_global$P.Value[match(pwy, sig_global$settings)],
        global_FDR = sig_global$adj.P.Val[match(pwy, sig_global$settings)],
        weight_boxplot_P = wp,
        LR_count_boxplot_P = cp
      )
    }
    safe_write_csv(bind_rows(master_box), file.path(plot_dirs$tables, "05_significant_pathway_boxplot_summary.csv"))
  }
}

# ------------------------------------------------------------
# 05. TARGET PATHWAYS TO PLOT FROM PER-SAMPLE OBJECTS
# ------------------------------------------------------------

all_pathways <- unique(unlist(lapply(obj_list, function(x) {
  if (!is.null(x@netP$pathways)) as.character(x@netP$pathways) else character(0)
})))

path_priority <- unique(c(TARGET_PATHWAYS, all_pathways))
path_priority <- path_priority[nzchar(path_priority)]
pathways_for_standard_plots <- head(path_priority, min(MAX_PATHWAYS_FOR_PER_SAMPLE, length(path_priority)))

# ------------------------------------------------------------
# 06. STANDARD CELLCHAT PER-SAMPLE NETWORK PLOTS
# ------------------------------------------------------------

if (MAKE_PER_SAMPLE_PLOTS) {
  samples_for_plots <- names(obj_list)
  if (is.finite(MAX_SAMPLES_FOR_FULL_PER_SAMPLE_PLOTS)) samples_for_plots <- head(samples_for_plots, MAX_SAMPLES_FOR_FULL_PER_SAMPLE_PLOTS)

  for (sid in samples_for_plots) {
    cc <- obj_list[[sid]]
    sdir <- file.path(plot_dirs$networks, sanitize(sid))
    pdir <- file.path(plot_dirs$cellchat_paths, sanitize(sid))
    dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

    # Overall interaction count and weight circles.
    if (!is.null(cc@net$count)) {
      safe_cellchat_plot(file.path(sdir, "01_circle_count.pdf"), function() {
        netVisual_circle(cc@net$count, vertex.weight = as.numeric(table(cc@idents)), weight.scale = TRUE,
                         label.edge = FALSE, title.name = paste0(sid, " - interaction count"))
      })
    }
    if (!is.null(cc@net$weight)) {
      safe_cellchat_plot(file.path(sdir, "02_circle_weight.pdf"), function() {
        netVisual_circle(cc@net$weight, vertex.weight = as.numeric(table(cc@idents)), weight.scale = TRUE,
                         label.edge = FALSE, title.name = paste0(sid, " - interaction weight"))
      })
    }

    # Standard CellChat bubble + heatmap for the full signaling network.
    safe_cellchat_plot(file.path(sdir, "03_bubble_all_interactions.pdf"), function() {
      netVisual_bubble(cc, remove.isolate = FALSE, title.name = paste0(sid, " - all interactions"))
    }, 10, 8)

    safe_cellchat_plot(file.path(sdir, "04_heatmap_all_signaling.pdf"), function() {
      netVisual_heatmap(cc, measure = "count", color.heatmap = "Reds", title.name = paste0(sid, " - signaling count"))
    }, 10, 10)

    safe_cellchat_plot(file.path(sdir, "05_heatmap_all_signaling_weight.pdf"), function() {
      netVisual_heatmap(cc, measure = "weight", color.heatmap = "Reds", title.name = paste0(sid, " - signaling weight"))
    }, 10, 10)

    # Signaling-role plots where available.
    available_pw_role <- if (!is.null(cc@netP$pathways)) as.character(cc@netP$pathways) else character(0)
    cc_role <- tryCatch(netAnalysis_computeCentrality(cc, slot.name = "netP"), error = function(e) NULL)

    if (!is.null(cc_role) && length(available_pw_role) > 0) {
      role_pw <- head(available_pw_role, min(20, length(available_pw_role)))
      safe_cellchat_plot(file.path(sdir, "06_signaling_role_sender_receiver.pdf"), function() {
        netAnalysis_signalingRole_network(cc_role, signaling = role_pw, width = 8, height = 2.5, font.size = 8)
      }, 10, 8)

      safe_cellchat_plot(file.path(sdir, "07_signaling_role_heatmap_incoming.pdf"), function() {
        netAnalysis_signalingRole_heatmap(cc_role, pattern = "incoming", signaling = role_pw, color.heatmap = "Reds")
      }, 10, 10)

      safe_cellchat_plot(file.path(sdir, "08_signaling_role_heatmap_outgoing.pdf"), function() {
        netAnalysis_signalingRole_heatmap(cc_role, pattern = "outgoing", signaling = role_pw, color.heatmap = "Reds")
      }, 10, 10)
    }

    # Pathway-specific standard CellChat plots.
    available_pw <- if (!is.null(cc@netP$pathways)) as.character(cc@netP$pathways) else character(0)
    pw_to_plot <- intersect(pathways_for_standard_plots, available_pw)

    for (pw in pw_to_plot) {
      z <- sanitize(pw)
      safe_cellchat_plot(file.path(pdir, paste0("01_circle_", z, ".pdf")), function() {
        netVisual_aggregate(cc, signaling = pw, layout = "circle", remove.isolate = FALSE)
      }, 8, 8)

      safe_cellchat_plot(file.path(pdir, paste0("02_hierarchy_", z, ".pdf")), function() {
        netVisual_aggregate(cc, signaling = pw, layout = "hierarchy", remove.isolate = FALSE)
      }, 9, 7)

      safe_cellchat_plot(file.path(pdir, paste0("03_heatmap_", z, ".pdf")), function() {
        netVisual_heatmap(cc, signaling = pw, color.heatmap = "Reds")
      }, 9, 8)

      safe_cellchat_plot(file.path(pdir, paste0("04_bubble_", z, ".pdf")), function() {
        netVisual_bubble(cc, signaling = pw, remove.isolate = FALSE, title.name = paste0(sid, " - ", pw))
      }, 10, 7)

      safe_cellchat_plot(file.path(pdir, paste0("05_contribution_", z, ".pdf")), function() {
        netAnalysis_contribution(cc, signaling = pw)
      }, 8, 6)
    }
  }
}

# ------------------------------------------------------------
# 07. GROUP-AVERAGED NETWORKS (DESCRIPTIVE, SAMPLE-AVERAGED)
# ------------------------------------------------------------

if (MAKE_GROUP_AVERAGE_NETWORKS && "AHI1_status" %in% colnames(coldata)) {
  sample_status <- coldata %>% select(Sample, AHI1_status) %>% distinct()
  sample_status <- sample_status %>% filter(Sample %in% names(obj_list))

  average_net <- function(objs, slot = c("count", "weight")) {
    slot <- match.arg(slot)
    mats <- lapply(objs, function(x) {
      m <- if (slot == "count") x@net$count else x@net$weight
      if (is.null(m)) return(NULL)
      as.matrix(m)
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (!length(mats)) return(NULL)
    common_sender <- Reduce(intersect, lapply(mats, rownames))
    common_receiver <- Reduce(intersect, lapply(mats, colnames))
    if (!length(common_sender) || !length(common_receiver)) return(NULL)
    mats2 <- lapply(mats, function(m) m[common_sender, common_receiver, drop = FALSE])
    Reduce(`+`, mats2) / length(mats2)
  }

  for (status in c("healthy", "risk")) {
    ids <- sample_status$Sample[sample_status$AHI1_status == status]
    ids <- intersect(ids, names(obj_list))
    if (!length(ids)) next
    gdir <- file.path(plot_dirs$group_networks, status)
    dir.create(gdir, recursive = TRUE, showWarnings = FALSE)

    for (metric in c("count", "weight")) {
      m <- average_net(obj_list[ids], metric)
      if (is.null(m)) next
      safe_cellchat_plot(file.path(gdir, paste0("01_average_network_", metric, ".pdf")), function() {
        netVisual_circle(m, vertex.weight = rowSums(m), weight.scale = TRUE, label.edge = FALSE,
                         title.name = paste0(status, " - sample-averaged ", metric))
      }, 8, 8)
      safe_write_csv(as.data.frame(as.table(m)) %>% rename(Sender = Var1, Receiver = Var2, value = Freq),
                     file.path(gdir, paste0("average_network_", metric, ".csv")))
    }
  }

  # Healthy-risks difference using sample-averaged networks.
  for (metric in c("count", "weight")) {
    mh <- average_net(obj_list[intersect(sample_status$Sample[sample_status$AHI1_status == "healthy"], names(obj_list))], metric)
    mr <- average_net(obj_list[intersect(sample_status$Sample[sample_status$AHI1_status == "risk"], names(obj_list))], metric)
    if (is.null(mh) || is.null(mr)) next
    common_r <- intersect(rownames(mh), rownames(mr)); common_c <- intersect(colnames(mh), colnames(mr))
    if (!length(common_r) || !length(common_c)) next
    d <- mr[common_r, common_c, drop = FALSE] - mh[common_r, common_c, drop = FALSE]
    safe_cellchat_plot(file.path(plot_dirs$group_networks, paste0("02_risk_minus_healthy_", metric, ".pdf")), function() {
      netVisual_circle(d, vertex.weight = rowSums(abs(d)), weight.scale = TRUE, label.edge = FALSE,
                       title.name = paste0("Risk - healthy average ", metric))
    }, 8, 8)
    safe_write_csv(as.data.frame(as.table(d)) %>% rename(Sender = Var1, Receiver = Var2, difference = Freq),
                   file.path(plot_dirs$group_networks, paste0("risk_minus_healthy_", metric, ".csv")))
  }
}

# ------------------------------------------------------------
# 08. TARGETED H1: GLUTAMATE AND ExN/InN -> OPC
# ------------------------------------------------------------

if (MAKE_TARGETED_PLOTS && "AHI1_status" %in% colnames(coldata)) {
  status <- coldata %>% select(Sample, AHI1_status) %>% distinct()

  h1_rows <- list()
  per_sample_h1 <- list()

  for (sid in intersect(names(obj_list), status$Sample)) {
    cc <- obj_list[[sid]]
    if (is.null(cc@net$prob)) next
    prob <- cc@net$prob
    dn <- dimnames(prob)
    if (length(dn) != 3) next
    senders <- intersect(H1_SENDERS, dn[[1]])
    if (!length(senders) || !(H1_RECEIVER %in% dn[[2]])) next
    ligands <- dn[[3]]

    glut_lrs <- tryCatch({
      db <- CellChatDB.human$interaction
      db$interaction_name[db$pathway_name == "Glutamate"]
    }, error = function(e) character(0))
    glut_lrs <- intersect(glut_lrs, ligands)
    if (!length(glut_lrs)) next

    sub <- prob[senders, H1_RECEIVER, glut_lrs, drop = FALSE]
    score <- sum(sub, na.rm = TRUE)
    h1_rows[[length(h1_rows) + 1L]] <- data.frame(
      Sample = sid,
      H1_glutamate_score = score
    )

    for (sender in senders) {
      ssub <- prob[sender, H1_RECEIVER, glut_lrs, drop = FALSE]
      per_sample_h1[[length(per_sample_h1) + 1L]] <- data.frame(
        Sample = sid,
        Sender = sender,
        Receiver = H1_RECEIVER,
        Pathway = "Glutamate",
        weight = sum(ssub, na.rm = TRUE)
      )
    }
  }

  h1 <- bind_rows(h1_rows) %>% left_join(status, by = "Sample")
  h1_sender <- bind_rows(per_sample_h1) %>% left_join(status, by = "Sample")

  if (nrow(h1) > 0) {
    safe_write_csv(h1, file.path(plot_dirs$targeted, "04_H1_ExN_InN_to_OPC_Glutamate_per_sample.csv"))
    p <- make_boxplot(h1, "H1_glutamate_score", "Aggregate ExN/InN -> OPC glutamatergic signaling", "Summed CellChat communication probability",
                      wilcox_safe(h1$H1_glutamate_score, h1$AHI1_status))
    if (!is.null(p)) safe_ggsave(file.path(plot_dirs$targeted, "05_H1_aggregate_Glutamate_ExN_InN_to_OPC.pdf"), p, 6, 5.5)
  }

  if (nrow(h1_sender) > 0) {
    p <- ggplot(h1_sender, aes(AHI1_status, weight, fill = AHI1_status)) +
      geom_boxplot(outlier.shape = NA, width = 0.55) +
      geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.75) +
      facet_wrap(~Sender, scales = "free_y") +
      scale_fill_manual(values = c(healthy = "#1b9e77", risk = "#d95f02")) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none") +
      labs(title = "Glutamate pathway: ExN/InN -> OPC", x = NULL, y = "Summed CellChat probability")
    safe_ggsave(file.path(plot_dirs$targeted, "06_Glutamate_ExN_vs_InN_to_OPC.pdf"), p, 8, 6)
    safe_write_csv(h1_sender, file.path(plot_dirs$targeted, "07_Glutamate_ExN_vs_InN_to_OPC_data.csv"))
  }

  # Targeted standard CellChat pathway plots for Glutamate / Ach.
  for (sid in intersect(names(obj_list), status$Sample)) {
    cc <- obj_list[[sid]]
    available_pw <- if (!is.null(cc@netP$pathways)) as.character(cc@netP$pathways) else character(0)
    for (pw in intersect(TARGET_PATHWAYS, available_pw)) {
      z <- sanitize(pw)
      safe_cellchat_plot(file.path(plot_dirs$targeted, paste0("08_", sanitize(sid), "_", z, "_ExN_InN_to_OPC_circle.pdf")), function() {
        src <- intersect(H1_SENDERS, rownames(cc@net$count))
        tgt <- if (H1_RECEIVER %in% colnames(cc@net$count)) H1_RECEIVER else NULL
        if (!length(src) || is.null(tgt)) stop("Targeted cell types not present in this CellChat object.")
        netVisual_aggregate(cc, signaling = pw, layout = "circle", sources.use = src, targets.use = tgt, remove.isolate = FALSE)
      }, 8, 8)
      safe_cellchat_plot(file.path(plot_dirs$targeted, paste0("09_", sanitize(sid), "_", z, "_bubble_ExN_InN_to_OPC.pdf")), function() {
        src <- intersect(H1_SENDERS, rownames(cc@net$count))
        tgt <- if (H1_RECEIVER %in% colnames(cc@net$count)) H1_RECEIVER else NULL
        if (!length(src) || is.null(tgt)) stop("Targeted cell types not present in this CellChat object.")
        netVisual_bubble(cc, sources.use = src, targets.use = tgt, signaling = pw, remove.isolate = FALSE,
                         title.name = paste0(sid, " - ", pw, " - ExN/InN -> OPC"))
      }, 9, 6)
    }
  }
}

# ------------------------------------------------------------
# 09. COPY EXISTING STEP-2 DIAGNOSTICS
# ------------------------------------------------------------

copy_names <- c(
  "LR_detection_glutamate_cholinergic.csv",
  "pathway_detection_glutamate_cholinergic.csv",
  "LR_detection_summary_before_filter.csv",
  "pathway_detection_summary_before_filter.csv",
  "LR_setting_detection_summary_before_filter.csv",
  "pathway_setting_detection_summary_before_filter.csv",
  "sample_metadata_used.csv"
)
for (f in copy_names) {
  src <- file.path(results_dir, f)
  if (file.exists(src)) file.copy(src, file.path(plot_dirs$tables, f), overwrite = TRUE)
}

# ------------------------------------------------------------
# 10. MANIFEST / SESSION INFO
# ------------------------------------------------------------

files_created <- list.files(plot_root, recursive = TRUE, full.names = FALSE)
safe_write_csv(data.frame(file = files_created), file.path(plot_dirs$manifest, "generated_files.csv"))

config_df <- data.frame(
  parameter = c(
    "BASE_OUT_DIR", "MODEL_NAME", "GROUP_DEFINITION", "RESTRICT_STAGE",
    "MIN_PROBABILITY", "MAX_TOP_LR", "MAX_TOP_PATHWAYS", "MAX_PATHWAYS_FOR_PER_SAMPLE",
    "MAKE_PER_SAMPLE_PLOTS", "MAKE_STANDARD_CELLCHAT_PLOTS", "MAKE_GROUP_AVERAGE_NETWORKS",
    "MAKE_TARGETED_PLOTS", "TARGET_PATHWAYS", "H1_SENDERS", "H1_RECEIVER"
  ),
  value = c(
    BASE_OUT_DIR, MODEL_NAME, GROUP_DEFINITION,
    ifelse(is.null(RESTRICT_STAGE), "NULL", RESTRICT_STAGE),
    MIN_PROBABILITY, MAX_TOP_LR, MAX_TOP_PATHWAYS, MAX_PATHWAYS_FOR_PER_SAMPLE,
    MAKE_PER_SAMPLE_PLOTS, MAKE_STANDARD_CELLCHAT_PLOTS, MAKE_GROUP_AVERAGE_NETWORKS,
    MAKE_TARGETED_PLOTS, paste(TARGET_PATHWAYS, collapse = ";"), paste(H1_SENDERS, collapse = ";"), H1_RECEIVER
  )
)
safe_write_csv(config_df, file.path(plot_dirs$manifest, "plot_configuration.csv"))

writeLines(capture.output(sessionInfo()), file.path(plot_dirs$manifest, "sessionInfo.txt"))

msg("============================================================")
msg("ALL CELLCHAT PLOTS COMPLETE")
msg("Plot directory: ", normalizePath(plot_root, mustWork = FALSE))
msg("Samples plotted: ", length(obj_list))
msg("============================================================")
