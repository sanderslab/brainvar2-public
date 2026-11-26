#!/usr/bin/env Rscript

# =============================================================================
# PHASE-LEVEL SINGLE-SCRIPT GRN PIPELINE
# =============================================================================
# Flow:
#   1) Phase-specific motif matching using phase-accessible peaks (all periods in phase)
#   2) Phase data retrieval by pooling period-level metacell matrices, then phase-level
#      batch correction across all period+batch groups
#   3) GRN inference (overlapped group lasso + post-selection inference) with unpenalized
#      period covariates
#
# Usage:
#   Rscript grn/run_phase_grn.R TARGET_GENE PHASE [REGION_SIZE] \
#     [FIMO_HIT_P_VALUE] [MIN_MEDIAN_OPEN_FRACTION] [PEAK_SCORE_QUANTILE] [N_RUNS] [P_VALUE_CUTOFF] \
#     [CELL_TYPE] [MAX_PHASE_PEAKS] [FIMO_SCAN_P] [N_CORES]
#
# Example:
#   Rscript grn/run_phase_grn.R BCL11B prenatal 1e5 1e-4 0.05 0.5 20 0.05 ExN 200 0.05 10

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(Matrix)
  library(biomaRt)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(progress)
  library(grpregOverlap)
})

safe_detect_cores <- function(default = 4L) {
  n <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(n) || is.na(n) || n < 1) return(as.integer(default))
  as.integer(n)
}

# -----------------------------------------------------------------------------
# Paths and source dependencies
# -----------------------------------------------------------------------------
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
INTERNAL_R_DIR <- file.path(SCRIPT_DIR, "R")
BASE_DIR <- SCRIPT_DIR
DATA_DIR <- Sys.getenv("GRN_DATA_DIR", unset = file.path(PROJECT_ROOT, "data"))
ATAC_DIR <- Sys.getenv("ATAC_DIR", unset = file.path(DATA_DIR, "atac"))
JASPAR_DIR <- Sys.getenv("MOTIF_DIR", unset = file.path(DATA_DIR, "JASPAR2024_CORE_pfms_meme"))
METACELL_DIR <- Sys.getenv(
  "METACELL_DIR",
  unset = file.path(PROJECT_ROOT, "explore", "explore_results", "metacells")
)
RESULTS_ROOT <- Sys.getenv("GRN_RESULTS_DIR", unset = file.path(PROJECT_ROOT, "results", "grn"))

ATAC_META_FILE <- Sys.getenv(
  "ATAC_METADATA_FILE",
  unset = file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-ATAC_P3-P15_meta_D240727.tsv")
)
MAIN_ARCHR_PROJ_FILE <- Sys.getenv(
  "ATAC_PROJECT_FILE",
  unset = file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-Save-ArchR-Project.rds")
)

PHASE_MOTIF_RESULTS_DIR <- file.path(RESULTS_ROOT, "motif")
PHASE_DATA_RETRIEVAL_DIR <- file.path(RESULTS_ROOT, "data")
PHASE_GRN_INFERENCE_DIR <- file.path(RESULTS_ROOT, "inference")
FIMO_FILES_DIR <- file.path(RESULTS_ROOT, "cache", "fasta")
FIMO_BG_DIR <- file.path(RESULTS_ROOT, "cache", "background")

PERIOD_MOTIF_SCRIPT <- file.path(INTERNAL_R_DIR, "motif_matching.R")
PERIOD_DATA_SCRIPT <- file.path(INTERNAL_R_DIR, "data_retrieval.R")
PERIOD_GRN_SCRIPT <- file.path(INTERNAL_R_DIR, "grn_inference.R")

for (f in c(PERIOD_MOTIF_SCRIPT, PERIOD_DATA_SCRIPT, PERIOD_GRN_SCRIPT)) {
  if (!file.exists(f)) stop("ERROR: Required pipeline script not found: ", f)
}

.SOURCING_PIPELINES <- TRUE
suppressMessages({
  capture.output({
    source(PERIOD_MOTIF_SCRIPT, local = TRUE)
    source(PERIOD_DATA_SCRIPT, local = TRUE)
    source(PERIOD_GRN_SCRIPT, local = TRUE)
  }, file = nullfile())
})
rm(.SOURCING_PIPELINES)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
normalize_peak_id <- function(x) {
  x <- gsub(":", "-", x)
  x <- gsub("\\.", "-", x)
  x
}

read_atac_meta <- function(path) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    return(as.data.frame(data.table::fread(path, sep = "auto", check.names = FALSE, data.table = FALSE)))
  }
  first_line <- readLines(path, n = 1, warn = FALSE)
  sep <- if (grepl(",", first_line, fixed = TRUE)) "," else "\t"
  if (sep == ",") {
    return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

cell_barcode_key <- function(x) {
  x <- as.character(x)
  x <- sub("^.*[_#]", "", x)
  x <- sub("-1$", "", x)
  x
}

extract_cell_sample <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("#", x), sub("#.*$", "", x), ifelse(grepl("_", x), sub("_.*$", "", x), NA_character_))
  out[!nzchar(out)] <- NA_character_
  out
}

canonical_cell_id <- function(x, sample_hint = NA_character_) {
  x <- as.character(x)
  sample <- extract_cell_sample(x)
  barcode <- cell_barcode_key(x)
  if (length(sample_hint) == 1) {
    sample_hint <- rep(sample_hint, length(x))
  } else {
    sample_hint <- as.character(sample_hint)
  }
  sample[is.na(sample) | !nzchar(sample)] <- sample_hint[is.na(sample) | !nzchar(sample)]
  ok <- !is.na(sample) & nzchar(sample) & nzchar(barcode)
  out <- rep(NA_character_, length(x))
  out[ok] <- paste0(sample[ok], "#", barcode[ok])
  out
}

detect_atac_meta_cell_column <- function(atac_meta) {
  candidates <- c("cellNames", "CellID", "cell_id", "Unnamed: 0", "X")
  hit <- candidates[candidates %in% colnames(atac_meta)]
  if (length(hit) > 0) {
    return(hit[1])
  }
  colnames(atac_meta)[1]
}

zscore_safe <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0 || !is.finite(m)) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

get_gene_info_from_coordinate_file <- function(target_gene) {
  coord_file <- Sys.getenv("GENE_COORDINATES_FILE", unset = file.path(PROJECT_ROOT, "config", "gene_coordinates.csv"))
  if (!file.exists(coord_file)) return(NULL)

  coords <- utils::read.csv(coord_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand")
  missing <- setdiff(required, colnames(coords))
  if (length(missing) > 0) {
    stop("Gene coordinate file is missing required columns: ", paste(missing, collapse = ", "))
  }

  hit <- coords[toupper(coords$hgnc_symbol) == toupper(target_gene), required, drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit[1, , drop = FALSE]
}

get_target_gene_info <- function(target_gene) {
  local_gene_info <- get_gene_info_from_coordinate_file(target_gene)
  if (!is.null(local_gene_info)) {
    cat("  Using local gene-coordinate file for", target_gene, "\n")
    return(local_gene_info)
  }

  gene_info <- tryCatch({
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    getBM(
      attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
      filters = "hgnc_symbol",
      values = target_gene,
      mart = mart
    )
  }, error = function(e) {
    cat("  Warning: Ensembl lookup failed for", target_gene, ":", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(gene_info) && nrow(gene_info) > 0) {
    return(gene_info[1, , drop = FALSE])
  }

  stop("Target gene not found in Ensembl or local gene-coordinate file: ", target_gene)
}

map_majorclass_to_broad <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  out[grepl("^ExN", x) | x %in% c("IT-Immature", "NonIT-Immature")] <- "ExN"
  out[grepl("^InN", x)] <- "InN"
  out[grepl("Glia", x)] <- "Glia"
  out[grepl("Progenitor", x)] <- "Progenitor"
  out[grepl("^NNC$", x)] <- "NNC"
  out
}

get_peak_ids_from_se <- function(se_obj) {
  rn <- rownames(se_obj)
  if (!is.null(rn) && length(rn) == nrow(se_obj)) {
    return(normalize_peak_id(rn))
  }
  rr <- tryCatch(rowRanges(se_obj), error = function(e) NULL)
  if (!is.null(rr) && length(rr) == nrow(se_obj)) {
    ids <- paste0(as.character(seqnames(rr)), "-", start(rr), "-", end(rr))
    return(normalize_peak_id(ids))
  }
  rep(NA_character_, nrow(se_obj))
}

get_valid_atac_batches_for_period <- function(period) {
  if (!file.exists(ATAC_META_FILE)) {
    stop("ERROR: ATAC metadata file not found: ", ATAC_META_FILE)
  }

  atac_meta <- read_atac_meta(ATAC_META_FILE)
  if (!("Sample" %in% colnames(atac_meta)) || !("Period" %in% colnames(atac_meta))) {
    stop("ERROR: ATAC metadata must contain 'Sample' and 'Period': ", ATAC_META_FILE)
  }

  period_num <- as.numeric(period)
  meta_period <- suppressWarnings(as.numeric(as.character(atac_meta$Period)))
  period_batches <- unique(as.character(atac_meta$Sample[!is.na(meta_period) & meta_period == period_num]))
  period_batches <- period_batches[nzchar(period_batches)]

  arrow_files <- list.files(
    ATAC_DIR,
    pattern = "ArrowFiles-.*\\.arrow$",
    full.names = FALSE
  )
  arrow_batches <- sub("\\.arrow$", "", sub("^.*ArrowFiles-", "", arrow_files))

  sort(intersect(period_batches, arrow_batches))
}

get_period_batches <- function(period) {
  pattern <- paste0("^mapping_p", period, "_.*\\.csv$")
  mapping_files <- list.files(METACELL_DIR, pattern = pattern, full.names = FALSE)
  if (length(mapping_files) == 0) return(character(0))

  map_batches <- gsub("\\.csv$", "", gsub(paste0("^mapping_p", period, "_"), "", mapping_files))
  valid_atac_batches <- get_valid_atac_batches_for_period(period)
  filtered <- sort(intersect(map_batches, valid_atac_batches))

  excluded <- setdiff(map_batches, filtered)
  if (length(excluded) > 0) {
    cat("Warning: Excluding mismatched/missing-Arrow batches for period", period, ":",
        paste(sort(excluded), collapse = ", "), "\n")
  }

  filtered
}

get_phase_periods <- function(phase) {
  phase <- tolower(phase)
  if (phase == "prenatal") return(3:6)
  if (phase == "postnatal") return(8:15)
  stop("ERROR: phase must be 'prenatal' or 'postnatal'")
}

get_phase_batch_table <- function(phase) {
  periods <- get_phase_periods(phase)
  out <- list()
  for (p in periods) {
    bs <- get_period_batches(p)
    if (length(bs) == 0) next
    out[[as.character(p)]] <- data.frame(period = p, batch = bs, stringsAsFactors = FALSE)
  }
  if (length(out) == 0) {
    return(data.frame(period = integer(0), batch = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, out)
}

get_archr_cells_for_phase <- function(phase_periods, batch_table, cell_type = NULL) {
  if (is.null(cell_type) || toupper(cell_type) == "ALL") {
    return(NULL)
  }
  atac_meta <- read_atac_meta(ATAC_META_FILE)
  cell_col <- detect_atac_meta_cell_column(atac_meta)
  if (!all(c("Sample", "Period", "rna_majorclass_v1") %in% colnames(atac_meta))) {
    stop("ATAC metadata missing required columns for cell-type filtering")
  }
  period_num <- suppressWarnings(as.numeric(as.character(atac_meta$Period)))
  broad <- map_majorclass_to_broad(atac_meta$rna_majorclass_v1)

  keep <- !is.na(period_num) & (period_num %in% as.numeric(phase_periods)) &
    (as.character(atac_meta$Sample) %in% batch_table$batch) &
    broad == cell_type

  unique(canonical_cell_id(
    x = as.character(atac_meta[[cell_col]][keep]),
    sample_hint = as.character(atac_meta$Sample[keep])
  ))
}

build_batch_metacell_groups <- function(batch_name, atac_cells, cell_type = NULL) {
  mc_assign <- load_metacell_assignments_filtered(batch_name, cell_type = cell_type)
  if (is.null(mc_assign) || length(mc_assign) == 0) {
    return(list())
  }

  map_cells <- unlist(mc_assign, use.names = FALSE)
  map_mcs <- rep(names(mc_assign), lengths(mc_assign))
  if (length(map_cells) == 0) {
    return(list())
  }

  mapping_df <- data.frame(
    key = canonical_cell_id(map_cells, sample_hint = batch_name),
    metacell = map_mcs,
    stringsAsFactors = FALSE
  )
  mapping_df <- mapping_df[!is.na(mapping_df$key) & nzchar(mapping_df$key), , drop = FALSE]
  if (nrow(mapping_df) == 0) {
    return(list())
  }

  atac_keys <- canonical_cell_id(atac_cells, sample_hint = batch_name)
  hit_idx <- match(atac_keys, mapping_df$key)
  valid <- !is.na(hit_idx)
  if (!any(valid)) {
    return(list())
  }

  split(atac_cells[valid], mapping_df$metacell[hit_idx[valid]])
}

safe_run <- function(expr, step_name) {
  tryCatch(expr, error = function(e) {
    stop("ERROR in ", step_name, ": ", e$message)
  })
}

motif_stage_complete <- function(stage_dir, target_gene) {
  req <- c(
    file.path(stage_dir, paste0(target_gene, "_motif_matrix.csv")),
    file.path(stage_dir, paste0(target_gene, "_metadata.rds")),
    file.path(stage_dir, paste0(target_gene, "_phase_peak_selection_stats.csv"))
  )
  all(file.exists(req))
}

motif_cache_matches <- function(stage_dir, target_gene, region_size,
                                fimo_hit_p_value_threshold,
                                min_median_open_fraction,
                                peak_score_quantile,
                                max_phase_peaks,
                                fimo_scan_p_value) {
  meta_file <- file.path(stage_dir, paste0(target_gene, "_metadata.rds"))
  if (!file.exists(meta_file)) return(FALSE)
  md <- tryCatch(readRDS(meta_file), error = function(e) NULL)
  if (is.null(md) || !is.list(md)) return(FALSE)

  checks <- c(
    isTRUE(all.equal(as.numeric(md$region_size), as.numeric(region_size), tolerance = 0)),
    isTRUE(all.equal(as.numeric(md$fimo_hit_p_value_threshold), as.numeric(fimo_hit_p_value_threshold), tolerance = 0)),
    isTRUE(all.equal(as.numeric(md$min_median_open_fraction), as.numeric(min_median_open_fraction), tolerance = 0)),
    isTRUE(all.equal(as.numeric(md$peak_score_quantile), as.numeric(peak_score_quantile), tolerance = 0)),
    isTRUE(all.equal(as.numeric(md$max_phase_peaks), as.numeric(max_phase_peaks), tolerance = 0)),
    isTRUE(all.equal(as.numeric(md$fimo_scan_p_value), as.numeric(fimo_scan_p_value), tolerance = 0))
  )
  all(checks)
}

data_stage_complete <- function(stage_dir, target_gene) {
  req <- c(
    file.path(stage_dir, paste0(target_gene, "_X_matrix.csv")),
    file.path(stage_dir, paste0(target_gene, "_Y_matrix.csv")),
    file.path(stage_dir, paste0(target_gene, "_Z_matrix.csv")),
    file.path(stage_dir, paste0(target_gene, "_motif_matrix.csv")),
    file.path(stage_dir, paste0(target_gene, "_batch_info.csv")),
    file.path(stage_dir, paste0(target_gene, "_data_retrieval_summary.txt"))
  )
  all(file.exists(req))
}

grn_stage_complete <- function(stage_dir, target_gene) {
  req <- c(
    file.path(stage_dir, paste0(target_gene, "_inference_summary.txt")),
    file.path(stage_dir, paste0(target_gene, "_combined_p_values.csv")),
    file.path(stage_dir, paste0(target_gene, "_significant_tfs.csv")),
    file.path(stage_dir, paste0(target_gene, "_significant_cres.csv")),
    file.path(stage_dir, paste0(target_gene, "_identified_tfs.csv")),
    file.path(stage_dir, paste0(target_gene, "_bootstrap_variable_details.csv")),
    file.path(stage_dir, paste0(target_gene, "_variable_score_summary.csv"))
  )
  all(file.exists(req))
}

grn_cache_matches <- function(stage_dir, target_gene, n_runs, p_value_cutoff) {
  summary_file <- file.path(stage_dir, paste0(target_gene, "_inference_summary.txt"))
  if (!file.exists(summary_file)) return(FALSE)
  lines <- tryCatch(readLines(summary_file, warn = FALSE), error = function(e) character(0))
  if (length(lines) == 0) return(FALSE)

  run_line <- grep("^Bootstrap Runs:", lines, value = TRUE)
  pval_line <- grep("^P-value Cutoff \\(BH-adjusted\\):", lines, value = TRUE)
  if (length(run_line) == 0 || length(pval_line) == 0) return(FALSE)

  cached_runs <- suppressWarnings(as.integer(sub("^Bootstrap Runs:\\s*", "", run_line[1])))
  cached_cutoff <- suppressWarnings(as.numeric(sub("^P-value Cutoff \\(BH-adjusted\\):\\s*", "", pval_line[1])))

  !is.na(cached_runs) &&
    !is.na(cached_cutoff) &&
    cached_runs == as.integer(n_runs) &&
    isTRUE(all.equal(cached_cutoff, as.numeric(p_value_cutoff), tolerance = 0))
}

load_cached_motif_results <- function(stage_dir, target_gene) {
  motif_matrix <- as.matrix(read.csv(
    file.path(stage_dir, paste0(target_gene, "_motif_matrix.csv")),
    row.names = 1,
    check.names = FALSE
  ))
  storage.mode(motif_matrix) <- "integer"
  metadata <- readRDS(file.path(stage_dir, paste0(target_gene, "_metadata.rds")) )
  peak_stats <- read.csv(
    file.path(stage_dir, paste0(target_gene, "_phase_peak_selection_stats.csv")),
    stringsAsFactors = FALSE
  )
  list(motif_matrix = motif_matrix, metadata = metadata, peak_stats = peak_stats, output_dir = stage_dir)
}

load_cached_data_results <- function(stage_dir, target_gene) {
  X_matrix <- as.matrix(read.csv(file.path(stage_dir, paste0(target_gene, "_X_matrix.csv")), row.names = 1, check.names = FALSE))
  Y_matrix <- as.matrix(read.csv(file.path(stage_dir, paste0(target_gene, "_Y_matrix.csv")), row.names = 1, check.names = FALSE))
  Z_matrix <- as.matrix(read.csv(file.path(stage_dir, paste0(target_gene, "_Z_matrix.csv")), row.names = 1, check.names = FALSE))
  motif_matrix <- as.matrix(read.csv(file.path(stage_dir, paste0(target_gene, "_motif_matrix.csv")), row.names = 1, check.names = FALSE))
  batch_info <- read.csv(file.path(stage_dir, paste0(target_gene, "_batch_info.csv")), stringsAsFactors = FALSE)
  list(X_matrix = X_matrix, Y_matrix = Y_matrix, Z_matrix = Z_matrix, motif_matrix = motif_matrix, batch_info = batch_info)
}

save_phase_grn_results <- function(results, target_gene, phase, output_dir, p_value_cutoff) {
  phase_output_dir <- file.path(output_dir, paste0("phase_", tolower(phase)), target_gene)
  if (!dir.exists(phase_output_dir)) {
    dir.create(phase_output_dir, recursive = TRUE)
  }

  sig_tfs_only_file <- file.path(phase_output_dir, paste0(target_gene, "_significant_tfs.csv"))
  if (length(results$significant_tfs_only) > 0) {
    score_idx <- match(results$significant_tfs_only, results$score_summary$Variable)
    sig_tfs_only_df <- data.frame(
      TF = results$significant_tfs_only,
      Mean_Score = results$score_summary$Mean_Coefficient[score_idx],
      Median_Score = results$score_summary$Median_Coefficient[score_idx],
      Selection_Frequency = results$score_summary$Selection_Frequency[score_idx],
      P_value = results$significant_tfs_only_pvals,
      stringsAsFactors = FALSE
    )
    sig_tfs_only_df <- sig_tfs_only_df[order(sig_tfs_only_df$P_value), ]
    write.csv(sig_tfs_only_df, sig_tfs_only_file, row.names = FALSE)
  } else {
    write.csv(data.frame(
      TF = character(0), Mean_Score = numeric(0), Median_Score = numeric(0),
      Selection_Frequency = numeric(0), P_value = numeric(0), stringsAsFactors = FALSE
    ), sig_tfs_only_file, row.names = FALSE)
  }

  sig_cres_file <- file.path(phase_output_dir, paste0(target_gene, "_significant_cres.csv"))
  if (length(results$significant_cres) > 0) {
    score_idx <- match(results$significant_cres, results$score_summary$Variable)
    sig_cres_df <- data.frame(
      CRE = results$significant_cres,
      Mean_Score = results$score_summary$Mean_Coefficient[score_idx],
      Median_Score = results$score_summary$Median_Coefficient[score_idx],
      Selection_Frequency = results$score_summary$Selection_Frequency[score_idx],
      P_value = results$significant_cres_pvals,
      stringsAsFactors = FALSE
    )
    sig_cres_df <- sig_cres_df[order(sig_cres_df$P_value), ]
    write.csv(sig_cres_df, sig_cres_file, row.names = FALSE)
  } else {
    write.csv(data.frame(
      CRE = character(0), Mean_Score = numeric(0), Median_Score = numeric(0),
      Selection_Frequency = numeric(0), P_value = numeric(0), stringsAsFactors = FALSE
    ), sig_cres_file, row.names = FALSE)
  }

  identified_tfs_file <- file.path(phase_output_dir, paste0(target_gene, "_identified_tfs.csv"))
  if (length(results$significant_tfs) > 0) {
    score_idx <- match(results$significant_tfs, results$score_summary$Variable)
    sig_tf_df <- data.frame(
      TF = results$significant_tfs,
      Mean_Score = results$score_summary$Mean_Coefficient[score_idx],
      Median_Score = results$score_summary$Median_Coefficient[score_idx],
      Selection_Frequency = results$score_summary$Selection_Frequency[score_idx],
      P_value = results$significant_tf_pvals,
      stringsAsFactors = FALSE
    )
    sig_tf_df <- sig_tf_df[order(sig_tf_df$P_value, na.last = TRUE), ]
    write.csv(sig_tf_df, identified_tfs_file, row.names = FALSE)
  } else {
    write.csv(data.frame(
      TF = character(0), Mean_Score = numeric(0), Median_Score = numeric(0),
      Selection_Frequency = numeric(0), P_value = numeric(0), stringsAsFactors = FALSE
    ), identified_tfs_file, row.names = FALSE)
  }

  pval_df <- data.frame(
    Variable = names(results$combined_p_values),
    P_value = results$combined_p_values,
    Adjusted_P_value = results$adjusted_p_values,
    stringsAsFactors = FALSE
  )
  write.csv(pval_df, file.path(phase_output_dir, paste0(target_gene, "_combined_p_values.csv")), row.names = FALSE)

  if (!is.null(results$bootstrap_diagnostics)) {
    diag_file <- file.path(phase_output_dir, paste0(target_gene, "_bootstrap_diagnostics.csv"))
    write.csv(results$bootstrap_diagnostics, diag_file, row.names = FALSE)
  }

  write.csv(
    results$bootstrap_details,
    file.path(phase_output_dir, paste0(target_gene, "_bootstrap_variable_details.csv")),
    row.names = FALSE
  )
  write.csv(
    results$score_summary,
    file.path(phase_output_dir, paste0(target_gene, "_variable_score_summary.csv")),
    row.names = FALSE
  )

  summary_file <- file.path(phase_output_dir, paste0(target_gene, "_inference_summary.txt"))
  summary_lines <- c(
    "Phase-Specific GRN Inference Pipeline Summary",
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Phase:", tolower(phase)),
    "Batch Correction (data retrieval): Yes (phase-level RPCA/LSI)",
    "Period Covariates (in model): Yes (unpenalized)",
    paste("Bootstrap Runs:", if (!is.null(results$p_val_matrix)) ncol(results$p_val_matrix) else NA_integer_),
    paste("P-value Cutoff (BH-adjusted):", p_value_cutoff),
    "Scores: mean/median coefficient across runs in which a variable was selected",
    "",
    "Results:",
    paste("  Significant TFs (p-value < cutoff):", ifelse(is.null(results$n_sig_tfs), length(results$significant_tfs_only), results$n_sig_tfs)),
    paste("  Significant CREs (p-value < cutoff):", ifelse(is.null(results$n_sig_cres), length(results$significant_cres), results$n_sig_cres)),
    paste("  Identified TFs (TF and binding CRE both significant):", ifelse(is.null(results$n_identified_tfs), length(results$significant_tfs), results$n_identified_tfs))
  )
  writeLines(summary_lines, summary_file)

  invisible(list(output_dir = phase_output_dir))
}

# -----------------------------------------------------------------------------
# Step 1: Phase-specific motif matching
# -----------------------------------------------------------------------------
run_phase_specific_motif_matching <- function(target_gene, phase, region_size = 1e5,
                                              fimo_hit_p_value_threshold = 1e-4,
                                              min_median_open_fraction = 0.05,
                                              peak_score_quantile = 0.5,
                                              cell_type = NULL,
                                              max_phase_peaks = 200,
                                              fimo_scan_p_value = 0.05,
                                              n_cores = 10,
                                              output_dir = PHASE_MOTIF_RESULTS_DIR) {
  cat("\n[1/3] Phase-specific motif matching\n")

  phase <- tolower(phase)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  phase_output_dir <- file.path(output_dir, paste0("phase_", phase), target_gene)
  if (!dir.exists(phase_output_dir)) dir.create(phase_output_dir, recursive = TRUE)

  period_vec <- get_phase_periods(phase)
  batch_table <- get_phase_batch_table(phase)
  if (nrow(batch_table) == 0) {
    stop("No valid ATAC batches found for phase ", phase, " after metadata+Arrow filtering")
  }

  cat("  Phase:", phase, "Periods:", paste(period_vec, collapse = ","), "\n")
  cat("  Valid period+batch pairs:", nrow(batch_table), "\n")
  cat("  Cell type:", ifelse(is.null(cell_type), "ALL", cell_type), "\n")

  if (!file.exists(MAIN_ARCHR_PROJ_FILE)) {
    stop("Main ArchR project not found: ", MAIN_ARCHR_PROJ_FILE)
  }
  proj <- readRDS(MAIN_ARCHR_PROJ_FILE)

  gene_info <- get_target_gene_info(target_gene)

  tss <- if (gene_info$strand == 1) gene_info$start_position else gene_info$end_position
  region <- GRanges(
    seqnames = paste0("chr", gene_info$chromosome_name),
    ranges = IRanges::IRanges(start = tss - region_size, end = tss + region_size)
  )

  peak_set <- getPeakSet(proj)
  if (is.null(peak_set) || length(peak_set) == 0) stop("No peak set in ArchR project")

  ov <- findOverlaps(peak_set, region)
  nearby_peaks <- peak_set[queryHits(ov)]
  if (length(nearby_peaks) == 0) stop("No peaks in target region for ", target_gene)

  nearby_coords <- paste(seqnames(nearby_peaks), start(nearby_peaks), end(nearby_peaks), sep = "-")
  near_norm <- normalize_peak_id(nearby_coords)
  cat("  Nearby peaks in region:", length(nearby_coords), "\n")

  allowed_cell_ids <- get_archr_cells_for_phase(period_vec, batch_table, cell_type)
  if (!is.null(allowed_cell_ids) && length(allowed_cell_ids) == 0) {
    cat("  Warning: No ATAC metadata keys found for selected cell type; using mapping-based matching only.\n")
    allowed_cell_ids <- NULL
  }

  batch_keys <- paste0("p", batch_table$period, "::", batch_table$batch)
  frac_open_mat <- matrix(
    NA_real_,
    nrow = length(nearby_coords),
    ncol = length(batch_keys),
    dimnames = list(nearby_coords, batch_keys)
  )
  mean_acc_mat <- matrix(
    NA_real_,
    nrow = length(nearby_coords),
    ncol = length(batch_keys),
    dimnames = list(nearby_coords, batch_keys)
  )

  cat("  Accessibility filter level: ATAC metacells across phase\n")

  for (i in seq_len(nrow(batch_table))) {
    p <- batch_table$period[i]
    b <- batch_table$batch[i]
    key <- batch_keys[i]

    proj_batch <- tryCatch(load_batch_archr_from_arrow(b), error = function(e) NULL)
    if (is.null(proj_batch)) next

    peak_se_b <- tryCatch(getMatrixFromProject(ArchRProj = proj_batch, useMatrix = "PeakMatrix"), error = function(e) NULL)
    if (is.null(peak_se_b)) next
    peak_mat_b <- assay(peak_se_b)
    row_norm_b <- get_peak_ids_from_se(peak_se_b)
    row_idx_b <- match(near_norm, row_norm_b)
    valid_rows <- !is.na(row_idx_b)
    if (!any(valid_rows)) next

    batch_cells <- colnames(peak_mat_b)
    if (!is.null(allowed_cell_ids)) {
      batch_cell_ids <- canonical_cell_id(batch_cells, sample_hint = b)
      keep_allowed <- batch_cell_ids %in% allowed_cell_ids
      if (any(keep_allowed)) {
        batch_cells <- batch_cells[keep_allowed]
      }
    }
    if (length(batch_cells) == 0) next

    mc_groups <- build_batch_metacell_groups(batch_name = b, atac_cells = batch_cells, cell_type = cell_type)
    if (length(mc_groups) == 0) next

    mc_groups <- mc_groups[lengths(mc_groups) > 0]
    if (length(mc_groups) == 0) next

    n_valid <- sum(valid_rows)
    n_mc <- length(mc_groups)
    open_bin <- matrix(0, nrow = n_valid, ncol = n_mc)
    mean_mc <- matrix(NA_real_, nrow = n_valid, ncol = n_mc)
    mat_peak <- peak_mat_b[row_idx_b[valid_rows], , drop = FALSE]

    j <- 0L
    for (mc_cells in mc_groups) {
      idx <- match(mc_cells, colnames(peak_mat_b))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) next
      j <- j + 1L
      sub_mat <- mat_peak[, idx, drop = FALSE]
      open_bin[, j] <- as.numeric(Matrix::rowSums(sub_mat > 0) > 0)
      mean_mc[, j] <- as.numeric(Matrix::rowMeans(sub_mat))
    }
    if (j == 0L) next

    frac_open_mat[valid_rows, key] <- rowMeans(open_bin[, seq_len(j), drop = FALSE], na.rm = TRUE)
    mean_acc_mat[valid_rows, key] <- rowMeans(mean_mc[, seq_len(j), drop = FALSE], na.rm = TRUE)
  }

  if (all(!is.finite(frac_open_mat))) {
    stop("No valid batch/cell values to compute phase-specific peak accessibility for phase ", phase)
  }

  batch_score_mat <- matrix(
    NA_real_,
    nrow = nrow(frac_open_mat),
    ncol = ncol(frac_open_mat),
    dimnames = dimnames(frac_open_mat)
  )

  for (b in colnames(frac_open_mat)) {
    frac_b <- frac_open_mat[, b]
    acc_b <- log1p(mean_acc_mat[, b])
    if (all(is.na(frac_b)) || all(is.na(acc_b))) next

    z_frac <- zscore_safe(frac_b)
    z_acc <- zscore_safe(acc_b)
    batch_score_mat[, b] <- 0.5 * z_frac + 0.5 * z_acc
  }

  phase_score <- rowMeans(batch_score_mat, na.rm = TRUE)
  median_frac_open <- apply(frac_open_mat, 1, median, na.rm = TRUE)
  mean_frac_open <- rowMeans(frac_open_mat, na.rm = TRUE)

  score_cutoff <- as.numeric(stats::quantile(phase_score, probs = peak_score_quantile, na.rm = TRUE, names = FALSE))
  keep <- is.finite(phase_score) & is.finite(median_frac_open) &
    (median_frac_open >= min_median_open_fraction) &
    (phase_score >= score_cutoff)

  if (!is.null(max_phase_peaks) && is.finite(max_phase_peaks) && max_phase_peaks > 0) {
    keep_idx <- which(keep)
    if (length(keep_idx) > max_phase_peaks) {
      ord <- order(phase_score[keep_idx], decreasing = TRUE, na.last = NA)
      keep_limited <- rep(FALSE, length(keep))
      keep_limited[keep_idx[ord[seq_len(max_phase_peaks)]]] <- TRUE
      keep <- keep_limited
    }
  }

  phase_specific_coords <- nearby_coords[keep]
  if (length(phase_specific_coords) == 0) {
    stop("No phase-specific peaks passed thresholds for ", target_gene, " in phase ", phase)
  }

  cat("  Phase-specific accessible peaks:", length(phase_specific_coords), "\n")
  cat("    Open-cell threshold (median across period+batches): >=", min_median_open_fraction, "\n")
  cat("    Batch-corrected peak score threshold (quantile): >=", peak_score_quantile,
      " (cutoff=", round(score_cutoff, 3), ")\n", sep = "")

  peak_stats <- data.frame(
    peak = nearby_coords,
    median_open_fraction = median_frac_open,
    mean_open_fraction = mean_frac_open,
    phase_score = phase_score,
    passed = keep,
    stringsAsFactors = FALSE
  )
  write.csv(peak_stats, file.path(phase_output_dir, paste0(target_gene, "_phase_peak_selection_stats.csv")), row.names = FALSE)

  if (!dir.exists(JASPAR_DIR)) stop("JASPAR dir not found: ", JASPAR_DIR)
  if (!dir.exists(FIMO_FILES_DIR)) dir.create(FIMO_FILES_DIR, recursive = TRUE)

  tf_collection <- extract_all_tf_names(JASPAR_DIR)
  if (length(tf_collection) == 0) stop("No TF motifs found in JASPAR directory")

  gr <- do.call(c, lapply(phase_specific_coords, function(x) {
    parts <- unlist(strsplit(x, "-"))
    GRanges(seqnames = parts[1], ranges = IRanges(as.numeric(parts[2]), as.numeric(parts[3])))
  }))
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  names(seqs) <- phase_specific_coords

  fasta_filepath <- file.path(FIMO_FILES_DIR, sprintf("%s_phase_%s_phase_specific_peaks.fa", target_gene, phase))
  writeXStringSet(seqs, filepath = fasta_filepath)

  bg_filename <- sprintf("%s_phase_%s_phase_specific_peaks.bg", target_gene, phase)
  prebuilt_bg <- file.path(FIMO_BG_DIR, bg_filename)
  bgfile_path <- if (file.exists(prebuilt_bg)) prebuilt_bg else file.path(FIMO_FILES_DIR, bg_filename)
  bgfile_built <- if (file.exists(prebuilt_bg)) {
    prebuilt_bg
  } else {
    build_markov_bgfile(fasta_filepath = fasta_filepath, bgfile_path = bgfile_path, order = 1L)
  }

  fimo_results_dir <- file.path(phase_output_dir, "fimo_results")
  if (!dir.exists(fimo_results_dir)) dir.create(fimo_results_dir, recursive = TRUE)

  total_tfs <- length(tf_collection)
  max_cores <- max(1L, min(as.integer(n_cores), safe_detect_cores(default = 4L), total_tfs))
  cat("  Running FIMO for", total_tfs, "TFs using", max_cores, "core(s)...\n")

  run_one_tf <- function(tf) {
    bvec <- tryCatch(
      run_fimo_pipeline_optimized(
        tf_folder = JASPAR_DIR,
        tf_name = tf,
        fimo_folder = fimo_results_dir,
        coords = phase_specific_coords,
        tg_name = target_gene,
        p_value_threshold = fimo_hit_p_value_threshold,
        fasta_filepath = fasta_filepath,
        scan_p_value_threshold = fimo_scan_p_value,
        bgfile_path = bgfile_built
      ),
      error = function(e) NULL
    )
    if (is.null(bvec)) {
      bvec <- rep(FALSE, length(phase_specific_coords))
    }
    list(tf = tf, bvec = bvec)
  }

  motif_library <- vector("list", total_tfs)
  names(motif_library) <- tf_collection

  if (max_cores > 1L && .Platform$OS.type != "windows") {
    chunk_size <- max(1L, max_cores * 2L)
    chunks <- split(seq_len(total_tfs), ceiling(seq_len(total_tfs) / chunk_size))
    done <- 0L
    for (chunk in chunks) {
      chunk_res <- parallel::mclapply(
        tf_collection[chunk],
        run_one_tf,
        mc.cores = max_cores,
        mc.preschedule = FALSE
      )
      for (j in seq_along(chunk)) {
        idx <- chunk[j]
        motif_library[[idx]] <- chunk_res[[j]]$bvec
      }
      done <- done + length(chunk)
      cat(sprintf("    FIMO progress: %d/%d TFs\n", done, total_tfs))
    }
  } else {
    pb <- progress_bar$new(total = total_tfs, format = "    [:bar] :current/:total (:percent) eta::eta", clear = FALSE, width = 60)
    for (i in seq_along(tf_collection)) {
      motif_library[[i]] <- run_one_tf(tf_collection[i])$bvec
      pb$tick()
    }
  }

  motif_matrix <- do.call(rbind, motif_library)
  colnames(motif_matrix) <- phase_specific_coords
  motif_matrix <- as.matrix(motif_matrix)
  storage.mode(motif_matrix) <- "integer"

  write.csv(motif_matrix, file.path(phase_output_dir, paste0(target_gene, "_motif_matrix.csv")))

  # Also write period-layout motif matrix copies for downstream period retrieval helper.
  phase_periods <- sort(unique(batch_table$period))
  for (p in phase_periods) {
    period_out_dir <- file.path(output_dir, paste0("period_", p), target_gene)
    if (!dir.exists(period_out_dir)) dir.create(period_out_dir, recursive = TRUE)
    write.csv(motif_matrix, file.path(period_out_dir, paste0(target_gene, "_motif_matrix.csv")))
  }

  metadata <- list(
    target_gene = target_gene,
    phase = phase,
    periods = phase_periods,
    region_size = region_size,
    fimo_hit_p_value_threshold = fimo_hit_p_value_threshold,
    min_median_open_fraction = min_median_open_fraction,
    peak_score_quantile = peak_score_quantile,
    max_phase_peaks = max_phase_peaks,
    fimo_scan_p_value = fimo_scan_p_value,
    score_cutoff = score_cutoff,
    tss = tss,
    region = as.character(region),
    n_nearby_peaks = length(nearby_coords),
    n_phase_specific_peaks = length(phase_specific_coords),
    n_tfs = nrow(motif_matrix),
    n_tfs_with_binding = sum(rowSums(motif_matrix) > 0),
    n_peaks_with_binding = sum(colSums(motif_matrix) > 0)
  )
  saveRDS(metadata, file.path(phase_output_dir, paste0(target_gene, "_metadata.rds")))

  cat("  Saved motif matrix:", file.path(phase_output_dir, paste0(target_gene, "_motif_matrix.csv")), "\n")
  cat("  Matrix shape:", nrow(motif_matrix), "TFs x", ncol(motif_matrix), "phase-specific peaks\n")

  list(motif_matrix = motif_matrix, metadata = metadata, peak_stats = peak_stats, output_dir = phase_output_dir)
}

# -----------------------------------------------------------------------------
# Step 2: Phase data retrieval and correction
# -----------------------------------------------------------------------------
run_phase_data_retrieval <- function(target_gene,
                                     phase,
                                     cell_type = "ALL",
                                     motif_results_root,
                                     output_dir = PHASE_DATA_RETRIEVAL_DIR,
                                     n_cores = 10) {
  cat("\n[2/3] Phase data retrieval\n")

  phase <- tolower(phase)
  phase_out <- file.path(output_dir, paste0("phase_", phase), target_gene)
  if (!dir.exists(phase_out)) dir.create(phase_out, recursive = TRUE)

  periods <- get_phase_periods(phase)
  raw_period_root <- file.path(phase_out, "period_raw")
  if (!dir.exists(raw_period_root)) dir.create(raw_period_root, recursive = TRUE)

  X_list <- list()
  Y_list <- list()
  Z_list <- list()
  info_list <- list()

  for (p in periods) {
    batches <- get_period_batches(p)
    if (length(batches) == 0) {
      cat("  Skipping period", p, "(no valid batch with ATAC Arrow + metacell mapping)\n")
      next
    }

    cat("  Retrieving period", p, "without intra-period correction...\n")
    per_res <- tryCatch(
      run_period_data_retrieval(
        target_gene = target_gene,
        period = p,
        output_dir = raw_period_root,
        apply_batch_correction = FALSE,
        cell_type = cell_type,
        motif_results_dir = motif_results_root,
        n_cores = n_cores
      ),
      error = function(e) {
        cat("    Warning: period", p, "failed:", e$message, "\n")
        NULL
      }
    )

    if (is.null(per_res)) next

    Xp <- per_res$X_matrix
    Yp <- per_res$Y_matrix
    Zp <- per_res$Z_matrix

    rid <- paste0("p", p, "::", rownames(Xp))
    rownames(Xp) <- rid
    rownames(Yp) <- rid
    rownames(Zp) <- rid

    X_list[[as.character(p)]] <- Xp
    Y_list[[as.character(p)]] <- Yp
    Z_list[[as.character(p)]] <- Zp

    bi <- data.frame(
      metacell_id = rid,
      batch = as.character(per_res$batch_info),
      period = as.integer(p),
      stringsAsFactors = FALSE
    )
    # Use period+batch key for cross-phase correction to avoid accidental label collision.
    bi$batch_phase <- paste0("p", bi$period, "::", bi$batch)
    info_list[[as.character(p)]] <- bi
  }

  if (length(X_list) == 0) {
    stop("No period produced valid X/Y/Z matrices for phase ", phase)
  }

  common_tfs <- Reduce(intersect, lapply(X_list, colnames))
  common_peaks <- Reduce(intersect, lapply(Z_list, colnames))
  if (length(common_tfs) == 0) stop("No common TF columns across periods in phase ", phase)
  if (length(common_peaks) == 0) stop("No common peak columns across periods in phase ", phase)

  X_combined <- do.call(rbind, lapply(X_list, function(x) x[, common_tfs, drop = FALSE]))
  Y_combined <- do.call(rbind, Y_list)
  Z_combined <- do.call(rbind, lapply(Z_list, function(z) z[, common_peaks, drop = FALSE]))
  Z_combined <- as.matrix(Z_combined)
  colnames(Z_combined) <- gsub("\\.", "-", colnames(Z_combined))

  batch_info_df <- do.call(rbind, info_list)
  batch_info_df <- batch_info_df[match(rownames(X_combined), batch_info_df$metacell_id), , drop = FALSE]
  if (any(is.na(batch_info_df$metacell_id))) {
    stop("Failed to align batch/period metadata with combined matrices")
  }

  x_cols_upper <- toupper(colnames(X_combined))
  target_upper <- toupper(target_gene)
  target_in_x <- which(x_cols_upper == target_upper)
  if (length(target_in_x) > 0) {
    X_combined <- X_combined[, -target_in_x, drop = FALSE]
    if (ncol(X_combined) == 0) stop("No X predictors after removing target gene")
  }

  cat("  Combined pre-correction matrices:\n")
  cat("    X:", dim(X_combined), "\n")
  cat("    Y:", dim(Y_combined), "\n")
  cat("    Z:", dim(Z_combined), "\n")

  n_unique_batches <- length(unique(batch_info_df$batch_phase))
  if (n_unique_batches <= 1) {
    cat("  Only one phase-batch found; skipping correction\n")
    X_corrected <- X_combined
    Y_corrected <- Y_combined
    Z_corrected <- Z_combined
  } else {
    cat("  Applying phase-level correction across", n_unique_batches, "period+batch groups (RNA RPCA, ATAC LSI)\n")

    temp_target_col <- "TARGETGENEFORRPCA"
    while (temp_target_col %in% colnames(X_combined)) {
      temp_target_col <- paste0(temp_target_col, "_X")
    }

    Y_for_correction <- Y_combined
    colnames(Y_for_correction) <- temp_target_col
    XY_combined <- cbind(X_combined, Y_for_correction)

    x_lookup <- colnames(X_combined)
    target_lookup <- temp_target_col
    if (anyDuplicated(colnames(XY_combined)) > 0) {
      colnames(XY_combined) <- make.unique(colnames(XY_combined), sep = "__dup")
      x_lookup <- colnames(XY_combined)[seq_len(ncol(X_combined))]
      target_lookup <- colnames(XY_combined)[ncol(XY_combined) + 1]
    }

    XY_corrected <- tryCatch(
      apply_rna_batch_correction(XY_combined, batch_info_df$batch_phase),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("max\\(nu, nv\\) must be strictly less than min\\(nrow\\(A\\), ncol\\(A\\)\\)", msg) ||
            grepl("Insufficient dimensions for RPCA", msg)) {
          cat("  Warning: RPCA skipped due to small sample size; using uncorrected RNA for this job.\n")
          return(XY_combined)
        }
        stop(e)
      }
    )

    xy_cols <- colnames(XY_corrected)
    if (is.null(xy_cols)) stop("Batch-corrected XY matrix has no column names")
    missing_x <- setdiff(x_lookup, xy_cols)
    if (length(missing_x) > 0) {
      stop("Missing TF columns after RNA correction: ", paste(head(missing_x, 5), collapse = ", "))
    }

    if (!(target_lookup %in% xy_cols)) {
      alt_target_lookup <- gsub("_", "-", target_lookup, fixed = TRUE)
      if (alt_target_lookup %in% xy_cols) {
        target_lookup <- alt_target_lookup
      } else {
        stop("Target feature column not found in corrected XY matrix")
      }
    }

    X_corrected <- XY_corrected[, x_lookup, drop = FALSE]
    colnames(X_corrected) <- colnames(X_combined)
    Y_corrected <- XY_corrected[, target_lookup, drop = FALSE]
    colnames(Y_corrected) <- target_gene

    Z_corrected <- apply_atac_batch_correction(Z_combined, batch_info_df$batch_phase)
  }

  motif_file <- file.path(motif_results_root, paste0("phase_", phase), target_gene, paste0(target_gene, "_motif_matrix.csv"))
  if (!file.exists(motif_file)) {
    stop("Phase motif matrix not found: ", motif_file)
  }
  motif_matrix <- as.matrix(read.csv(motif_file, row.names = 1, check.names = FALSE))

  write.csv(X_corrected, file.path(phase_out, paste0(target_gene, "_X_matrix.csv")))
  write.csv(Y_corrected, file.path(phase_out, paste0(target_gene, "_Y_matrix.csv")))
  write.csv(Z_corrected, file.path(phase_out, paste0(target_gene, "_Z_matrix.csv")))
  write.csv(motif_matrix, file.path(phase_out, paste0(target_gene, "_motif_matrix.csv")))
  write.csv(batch_info_df, file.path(phase_out, paste0(target_gene, "_batch_info.csv")), row.names = FALSE)

  summary_lines <- c(
    "Phase-Specific Data Retrieval Pipeline Summary",
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Phase:", phase),
    paste("Cell Type:", ifelse(is.null(cell_type), "ALL", cell_type)),
    paste("Periods Included:", paste(sort(unique(batch_info_df$period)), collapse = ", ")),
    paste("Unique period+batch groups:", n_unique_batches),
    "",
    "Matrix Dimensions:",
    paste("  X (TFs):", nrow(X_corrected), "metacells x", ncol(X_corrected), "TFs"),
    paste("  Y (Target):", nrow(Y_corrected), "metacells x", ncol(Y_corrected), "target"),
    paste("  Z (ATAC):", nrow(Z_corrected), "metacells x", ncol(Z_corrected), "peaks"),
    "",
    "Batch Correction:",
    "  RNA (X, Y): Seurat RPCA across all period+batch groups",
    "  ATAC (Z): LSI embedding correction across all period+batch groups",
    paste("Output Directory:", phase_out)
  )
  writeLines(summary_lines, file.path(phase_out, paste0(target_gene, "_data_retrieval_summary.txt")))

  list(
    X_matrix = X_corrected,
    Y_matrix = Y_corrected,
    Z_matrix = Z_corrected,
    motif_matrix = motif_matrix,
    batch_info = batch_info_df,
    output_dir = phase_out
  )
}

# -----------------------------------------------------------------------------
# Step 3: Phase GRN inference with unpenalized period covariates
# -----------------------------------------------------------------------------
run_phase_grn_inference <- function(target_gene,
                                    phase,
                                    data_results,
                                    n_runs = 100,
                                    p_value_cutoff = 0.05,
                                    output_dir = PHASE_GRN_INFERENCE_DIR) {
  cat("\n[3/3] Phase GRN inference\n")

  aligned <- verify_matrices(list(
    X_matrix = data_results$X_matrix,
    Y_matrix = data_results$Y_matrix,
    Z_matrix = data_results$Z_matrix,
    motif_matrix = data_results$motif_matrix,
    batch_info = data_results$batch_info
  ))

  x_cols_upper <- toupper(colnames(aligned$X_matrix))
  target_upper <- toupper(target_gene)
  target_in_x_idx <- which(x_cols_upper == target_upper)
  if (length(target_in_x_idx) > 0) {
    aligned$X_matrix <- aligned$X_matrix[, -target_in_x_idx, drop = FALSE]
    cat("  Removed target gene from X covariates before GRN fit\n")
  }
  if (ncol(aligned$X_matrix) == 0) {
    stop("ERROR: X matrix has no predictors after removing target gene")
  }

  period_covariates <- NULL
  if (!is.null(aligned$batch_info) && ("period" %in% colnames(aligned$batch_info))) {
    period_labels <- as.factor(aligned$batch_info$period)
    if (nlevels(period_labels) > 1) {
      mm <- model.matrix(~ period_labels)
      if (ncol(mm) > 1) {
        period_covariates <- mm[, -1, drop = FALSE]
        colnames(period_covariates) <- paste0("period_", levels(period_labels)[-1])
      }
    }
  }

  phase_out <- file.path(output_dir, paste0("phase_", tolower(phase)), target_gene)
  if (!dir.exists(phase_out)) dir.create(phase_out, recursive = TRUE)
  progress_file <- file.path(phase_out, paste0(target_gene, "_bootstrap_progress.txt"))

  grn_results <- run_grn_inference(
    aligned$X_matrix,
    aligned$Y_matrix,
    aligned$Z_matrix,
    aligned$motif_matrix,
    batch_covariates = period_covariates,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff,
    progress_file = progress_file,
    progress_every = 5
  )

  save_phase_grn_results(
    results = grn_results,
    target_gene = target_gene,
    phase = phase,
    output_dir = output_dir,
    p_value_cutoff = p_value_cutoff
  )

  grn_results
}

# -----------------------------------------------------------------------------
# Main single-script pipeline
# -----------------------------------------------------------------------------
run_phase_grn_pipeline <- function(target_gene,
                                              phase,
                                              region_size = 1e5,
                                              fimo_hit_p_value_threshold = 1e-4,
                                              min_median_open_fraction = 0.05,
                                              peak_score_quantile = 0.5,
                                              n_runs = 100,
                                              p_value_cutoff = 0.05,
                                              cell_type = "ALL",
                                              max_phase_peaks = 200,
                                              fimo_scan_p_value = 0.05,
                                              n_cores = 10) {
  phase <- tolower(phase)

  cat("\n", paste(rep("=", 88), collapse = ""), "\n", sep = "")
  cat("PHASE-LEVEL SINGLE-SCRIPT GRN PIPELINE\n")
  cat(paste(rep("=", 88), collapse = ""), "\n")
  cat("Target gene:", target_gene, "\n")
  cat("Phase:", phase, "\n")
  cat("Region size:", region_size, "bp\n")
  cat("FIMO hit p-value threshold:", fimo_hit_p_value_threshold, "\n")
  cat("FIMO scan p-value threshold:", fimo_scan_p_value, "\n")
  cat("Minimum median open-cell fraction:", min_median_open_fraction, "\n")
  cat("Peak score quantile:", peak_score_quantile, "\n")
  cat("Maximum retained phase peaks:", max_phase_peaks, "\n")
  cat("Parallel cores:", n_cores, "\n")
  cat("GRN bootstrap runs:", n_runs, "\n")
  cat("GRN p-value cutoff:", p_value_cutoff, "\n")
  cat("Batch correction (fixed): TRUE (phase-level)\n")
  cat("Regression covariates (fixed): period indicators, unpenalized\n")
  cat("Cell type:", cell_type, "\n")

  cell_type_suffix <- if (is.null(cell_type) || toupper(cell_type) == "ALL") "all_cells" else paste0("celltype_", cell_type)

  motif_base <- file.path(PHASE_MOTIF_RESULTS_DIR, cell_type_suffix)
  data_base <- file.path(PHASE_DATA_RETRIEVAL_DIR, cell_type_suffix)
  grn_base <- file.path(PHASE_GRN_INFERENCE_DIR, cell_type_suffix)

  motif_stage_dir <- file.path(motif_base, paste0("phase_", phase), target_gene)
  data_stage_dir <- file.path(data_base, paste0("phase_", phase), target_gene)
  grn_stage_dir <- file.path(grn_base, paste0("phase_", phase), target_gene)

  addArchRThreads(threads = min(4, max(1, as.integer(n_cores))))
  addArchRGenome("hg38")

  motif_ready <- motif_stage_complete(motif_stage_dir, target_gene) &&
    motif_cache_matches(
      stage_dir = motif_stage_dir,
      target_gene = target_gene,
      region_size = region_size,
      fimo_hit_p_value_threshold = fimo_hit_p_value_threshold,
      min_median_open_fraction = min_median_open_fraction,
      peak_score_quantile = peak_score_quantile,
      max_phase_peaks = max_phase_peaks,
      fimo_scan_p_value = fimo_scan_p_value
    )
  data_ready <- data_stage_complete(data_stage_dir, target_gene)
  grn_ready <- grn_stage_complete(grn_stage_dir, target_gene) &&
    grn_cache_matches(grn_stage_dir, target_gene, n_runs = n_runs, p_value_cutoff = p_value_cutoff)

  if (grn_ready && data_ready && motif_ready) {
    cat("\n[Cache] Completed phase GRN results already exist for this target/cell-type/phase.\n")
    return(list(motif_results = NULL, data_results = NULL, grn_results = list(output_dir = grn_stage_dir, cached = TRUE)))
  }

  if (motif_ready) {
    cat("\n[Cache] Using existing phase motif results from:\n  ", motif_stage_dir, "\n", sep = "")
    motif_results <- load_cached_motif_results(motif_stage_dir, target_gene)
  } else {
    motif_results <- safe_run(
      run_phase_specific_motif_matching(
        target_gene = target_gene,
        phase = phase,
        region_size = region_size,
        fimo_hit_p_value_threshold = fimo_hit_p_value_threshold,
        min_median_open_fraction = min_median_open_fraction,
        peak_score_quantile = peak_score_quantile,
        cell_type = cell_type,
        max_phase_peaks = max_phase_peaks,
        fimo_scan_p_value = fimo_scan_p_value,
        n_cores = n_cores,
        output_dir = motif_base
      ),
      "phase-specific motif matching"
    )
  }

  if (data_ready) {
    cat("\n[Cache] Using existing phase data-retrieval results from:\n  ", data_stage_dir, "\n", sep = "")
    data_results <- load_cached_data_results(data_stage_dir, target_gene)
  } else {
    data_results <- safe_run(
      run_phase_data_retrieval(
        target_gene = target_gene,
        phase = phase,
        cell_type = cell_type,
        motif_results_root = motif_base,
        output_dir = data_base,
        n_cores = n_cores
      ),
      "phase data retrieval"
    )
  }

  if (grn_ready) {
    cat("\n[Cache] Using existing phase GRN-inference results from:\n  ", grn_stage_dir, "\n", sep = "")
    grn_results <- list(output_dir = grn_stage_dir, cached = TRUE)
  } else {
    grn_results <- safe_run(
      run_phase_grn_inference(
        target_gene = target_gene,
        phase = phase,
        data_results = data_results,
        n_runs = n_runs,
        p_value_cutoff = p_value_cutoff,
        output_dir = grn_base
      ),
      "phase GRN inference"
    )
  }

  cat("\n", paste(rep("=", 88), collapse = ""), "\n", sep = "")
  cat("PIPELINE COMPLETE\n")
  cat(paste(rep("=", 88), collapse = ""), "\n")
  cat("Motif results:", motif_stage_dir, "\n")
  cat("Data retrieval results:", data_stage_dir, "\n")
  cat("GRN inference results:", grn_stage_dir, "\n")

  list(motif_results = motif_results, data_results = data_results, grn_results = grn_results)
}

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
if (length(commandArgs(trailingOnly = TRUE)) > 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2) {
    stop("Usage: Rscript grn/run_phase_grn.R TARGET_GENE PHASE [REGION_SIZE] [FIMO_HIT_P_VALUE] [MIN_MEDIAN_OPEN_FRACTION] [PEAK_SCORE_QUANTILE] [N_RUNS] [P_VALUE_CUTOFF] [CELL_TYPE] [MAX_PHASE_PEAKS] [FIMO_SCAN_P] [N_CORES]")
  }

  target_gene <- args[1]
  phase <- tolower(args[2])
  region_size <- if (length(args) > 2) as.numeric(args[3]) else 1e5
  fimo_hit_p_value_threshold <- if (length(args) > 3) as.numeric(args[4]) else 1e-4
  min_median_open_fraction <- if (length(args) > 4) as.numeric(args[5]) else 0.05
  peak_score_quantile <- if (length(args) > 5) as.numeric(args[6]) else 0.5
  n_runs <- if (length(args) > 6) as.integer(args[7]) else 100
  p_value_cutoff <- if (length(args) > 7) as.numeric(args[8]) else 0.05
  cell_type <- if (length(args) > 8) args[9] else "ALL"
  max_phase_peaks <- if (length(args) > 9) as.integer(args[10]) else 200
  fimo_scan_p_value <- if (length(args) > 10) as.numeric(args[11]) else 0.05
  n_cores <- if (length(args) > 11) as.integer(args[12]) else 10

  if (!(phase %in% c("prenatal", "postnatal"))) stop("ERROR: phase must be one of prenatal, postnatal")
  if (is.na(region_size) || region_size <= 0) stop("ERROR: region_size must be > 0")
  if (is.na(fimo_hit_p_value_threshold) || fimo_hit_p_value_threshold <= 0 || fimo_hit_p_value_threshold > 1) stop("ERROR: fimo_hit_p_value_threshold must be in (0,1]")
  if (is.na(min_median_open_fraction) || min_median_open_fraction < 0 || min_median_open_fraction > 1) stop("ERROR: min_median_open_fraction must be in [0,1]")
  if (is.na(peak_score_quantile) || peak_score_quantile < 0 || peak_score_quantile > 1) stop("ERROR: peak_score_quantile must be in [0,1]")
  if (!is.na(max_phase_peaks) && max_phase_peaks < 1) stop("ERROR: max_phase_peaks must be >= 1")
  if (is.na(fimo_scan_p_value) || fimo_scan_p_value <= 0 || fimo_scan_p_value > 1) stop("ERROR: fimo_scan_p_value must be in (0,1]")
  if (is.na(n_runs) || n_runs <= 0) stop("ERROR: n_runs must be a positive integer")
  if (is.na(n_cores) || n_cores < 1) stop("ERROR: n_cores must be a positive integer")
  if (!(toupper(cell_type) %in% c("ALL", "EXN", "INN", "GLIA", "PROGENITOR", "NNC"))) {
    stop("ERROR: CELL_TYPE must be one of ALL, ExN, InN, Glia, Progenitor, NNC")
  }

  if (toupper(cell_type) != "ALL") {
    cell_type <- c("EXN" = "ExN", "INN" = "InN", "GLIA" = "Glia", "PROGENITOR" = "Progenitor", "NNC" = "NNC")[toupper(cell_type)]
  } else {
    cell_type <- "ALL"
  }

  run_phase_grn_pipeline(
    target_gene = target_gene,
    phase = phase,
    region_size = region_size,
    fimo_hit_p_value_threshold = fimo_hit_p_value_threshold,
    min_median_open_fraction = min_median_open_fraction,
    peak_score_quantile = peak_score_quantile,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff,
    cell_type = cell_type,
    max_phase_peaks = max_phase_peaks,
    fimo_scan_p_value = fimo_scan_p_value,
    n_cores = n_cores
  )
} else {
  cat("Usage: Rscript grn/run_phase_grn.R TARGET_GENE PHASE [REGION_SIZE] [FIMO_HIT_P_VALUE] [MIN_MEDIAN_OPEN_FRACTION] [PEAK_SCORE_QUANTILE] [N_RUNS] [P_VALUE_CUTOFF] [CELL_TYPE] [MAX_PHASE_PEAKS] [FIMO_SCAN_P] [N_CORES]\n")
}
