#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(IRanges)
  library(Matrix)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(biomaRt)
})

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1])), mustWork = TRUE))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

SCRIPT_DIR <- get_script_dir()
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
source(file.path(SCRIPT_DIR, "R", "grn_scoring.R"))
DATA_DIR <- Sys.getenv("GRN_DATA_DIR", unset = file.path(PROJECT_ROOT, "data"))
ATAC_DIR <- Sys.getenv("ATAC_DIR", unset = file.path(DATA_DIR, "atac"))
JASPAR_DIR <- Sys.getenv("MOTIF_DIR", unset = file.path(DATA_DIR, "JASPAR2024_CORE_pfms_meme"))
METACELL_DIR <- Sys.getenv(
  "METACELL_DIR",
  unset = file.path(PROJECT_ROOT, "explore", "explore_results", "metacells")
)
RESULTS_ROOT <- Sys.getenv("COLOC_RESULTS_DIR", unset = file.path(PROJECT_ROOT, "results", "coloc"))
LOG_ROOT <- file.path(RESULTS_ROOT, "logs")
TMP_ROOT <- file.path(RESULTS_ROOT, "tmp")
MAIN_ARCHR_PROJ_FILE <- Sys.getenv(
  "ATAC_PROJECT_FILE",
  unset = file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-Save-ArchR-Project.rds")
)
ATAC_META_FILE <- Sys.getenv(
  "ATAC_METADATA_FILE",
  unset = file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-ATAC_P3-P15_meta_D240727.tsv")
)

dir.create(RESULTS_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_ROOT, recursive = TRUE, showWarnings = FALSE)
dir.create(TMP_ROOT, recursive = TRUE, showWarnings = FALSE)

source_pipeline_helpers <- function(project_root) {
  helpers_env <- new.env(parent = globalenv())
  helpers_env$.SOURCING_PIPELINES <- TRUE
  helpers_env$PROJECT_ROOT <- project_root
  helpers_env$BASE_DIR <- file.path(project_root, "grn")
  helpers_env$INTERNAL_R_DIR <- file.path(project_root, "grn", "R")
  helpers_env$DATA_DIR <- DATA_DIR
  helpers_env$ATAC_DIR <- ATAC_DIR
  helpers_env$JASPAR_DIR <- JASPAR_DIR
  helpers_env$METACELL_DIR <- METACELL_DIR
  helpers_env$MAIN_ARCHR_PROJ_FILE <- MAIN_ARCHR_PROJ_FILE
  helpers_env$ATAC_META_FILE <- ATAC_META_FILE
  source(file.path(helpers_env$INTERNAL_R_DIR, "motif_matching.R"), local = helpers_env)
  source(file.path(helpers_env$INTERNAL_R_DIR, "data_retrieval.R"), local = helpers_env)
  source(file.path(helpers_env$INTERNAL_R_DIR, "grn_inference.R"), local = helpers_env)
  helpers_env
}

helpers <- source_pipeline_helpers(PROJECT_ROOT)

normalize_peak_id <- function(x) {
  x <- gsub(":", "-", x)
  gsub("\\.", "-", x)
}

detect_atac_meta_cell_column <- function(atac_meta) {
  candidates <- c("cellNames", "CellID", "cell_id", "Unnamed: 0", "X")
  hit <- candidates[candidates %in% colnames(atac_meta)]
  if (length(hit) > 0) {
    return(hit[1])
  }
  colnames(atac_meta)[1]
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

zscore_safe <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0 || !is.finite(m)) {
    return(rep(0, length(x)))
  }
  (x - m) / s
}

scale_matrix_columns <- function(mat, log1p_nonnegative = FALSE) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  if (log1p_nonnegative) {
    mat <- log1p(pmax(mat, 0))
  }
  out <- apply(mat, 2, zscore_safe)
  if (is.null(dim(out))) {
    out <- matrix(out, ncol = 1)
  }
  out <- as.matrix(out)
  rownames(out) <- rownames(mat)
  colnames(out) <- colnames(mat)
  out
}

safe_make_names <- function(x) {
  make.names(as.character(x), unique = TRUE)
}

coerce_chr <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^chr", x, ignore.case = TRUE), x, paste0("chr", x))
}

read_loci <- function(path) {
  if (!file.exists(path)) {
    stop("Loci file not found: ", path)
  }
  loci_df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("rsid", "chr", "pos") %in% colnames(loci_df))) {
    stop("Loci CSV must contain columns: rsid, chr, pos")
  }
  loci_df$chr <- coerce_chr(loci_df$chr)
  loci_df$pos <- as.integer(loci_df$pos)
  loci_df <- loci_df[!is.na(loci_df$pos) & nzchar(loci_df$chr), , drop = FALSE]
  if (nrow(loci_df) == 0) {
    stop("No valid loci rows found in: ", path)
  }
  loci_df
}

get_phase_periods <- function(phase) {
  phase <- tolower(phase)
  if (phase == "prenatal") {
    return(3:6)
  }
  if (phase == "postnatal") {
    return(8:15)
  }
  stop("Phase must be prenatal or postnatal")
}

get_valid_atac_batches_for_period <- function(period) {
  atac_meta <- helpers$read_atac_meta(ATAC_META_FILE)
  if (!all(c("Sample", "Period") %in% colnames(atac_meta))) {
    stop("ATAC metadata must contain Sample and Period columns")
  }
  meta_period <- suppressWarnings(as.numeric(as.character(atac_meta$Period)))
  period_batches <- unique(as.character(atac_meta$Sample[!is.na(meta_period) & meta_period == as.numeric(period)]))
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
  mapping_dir <- METACELL_DIR
  pattern <- paste0("^mapping_p", period, "_.*\\.csv$")
  mapping_files <- list.files(mapping_dir, pattern = pattern, full.names = FALSE)
  if (length(mapping_files) == 0) {
    return(character(0))
  }
  map_batches <- gsub("\\.csv$", "", gsub(paste0("^mapping_p", period, "_"), "", mapping_files))
  valid_atac_batches <- get_valid_atac_batches_for_period(period)
  filtered <- sort(intersect(map_batches, valid_atac_batches))
  excluded <- setdiff(map_batches, filtered)
  if (length(excluded) > 0) {
    cat("Warning: excluding mismatched/missing-Arrow batches for period", period, ":",
        paste(excluded, collapse = ", "), "\n")
  }
  filtered
}

get_phase_batch_table <- function(phase) {
  periods <- get_phase_periods(phase)
  batch_rows <- list()
  for (p in periods) {
    batches <- get_period_batches(p)
    if (length(batches) == 0) {
      next
    }
    batch_rows[[as.character(p)]] <- data.frame(
      period = as.integer(p),
      batch = batches,
      stringsAsFactors = FALSE
    )
  }
  if (length(batch_rows) == 0) {
    return(data.frame(period = integer(0), batch = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, batch_rows)
}

get_archr_cells_for_phase <- function(phase_periods, batch_table, cell_type = NULL) {
  if (is.null(cell_type) || toupper(cell_type) == "ALL") {
    return(NULL)
  }
  atac_meta <- helpers$read_atac_meta(ATAC_META_FILE)
  cell_col <- detect_atac_meta_cell_column(atac_meta)
  required_cols <- c("Sample", "Period", "rna_majorclass_v1")
  if (!all(required_cols %in% colnames(atac_meta))) {
    stop("ATAC metadata missing required columns for cell-type filtering")
  }
  period_num <- suppressWarnings(as.numeric(as.character(atac_meta$Period)))
  broad <- helpers$map_majorclass_to_broad(atac_meta$rna_majorclass_v1)
  keep <- !is.na(period_num) &
    (period_num %in% as.numeric(phase_periods)) &
    (as.character(atac_meta$Sample) %in% batch_table$batch) &
    broad == cell_type

  unique(helpers$canonical_cell_id(
    as.character(atac_meta[[cell_col]][keep]),
    sample_hint = as.character(atac_meta$Sample[keep])
  ))
}

build_batch_metacell_groups <- function(batch_name, atac_cells, cell_type = NULL) {
  mc_assign <- helpers$load_metacell_assignments_filtered(batch_name, cell_type = cell_type)
  if (is.null(mc_assign) || length(mc_assign) == 0) {
    return(list())
  }

  map_cells <- unlist(mc_assign, use.names = FALSE)
  map_mcs <- rep(names(mc_assign), lengths(mc_assign))
  if (length(map_cells) == 0) {
    return(list())
  }

  mapping_df <- data.frame(
    key = helpers$canonical_cell_id(map_cells, sample_hint = batch_name),
    metacell = map_mcs,
    stringsAsFactors = FALSE
  )
  mapping_df <- mapping_df[!is.na(mapping_df$key) & nzchar(mapping_df$key), , drop = FALSE]
  if (nrow(mapping_df) == 0) {
    return(list())
  }

  atac_keys <- helpers$canonical_cell_id(atac_cells, sample_hint = batch_name)
  hit_idx <- match(atac_keys, mapping_df$key)
  valid <- !is.na(hit_idx)
  if (!any(valid)) {
    return(list())
  }
  split(atac_cells[valid], mapping_df$metacell[hit_idx[valid]])
}

fit_lm_term <- function(response, predictor, predictor_name, covariates = NULL) {
  response <- as.numeric(response)
  predictor <- as.numeric(predictor)

  if (!is.finite(stats::sd(response, na.rm = TRUE)) || stats::sd(response, na.rm = TRUE) == 0) {
    return(list(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = length(response)))
  }
  if (!is.finite(stats::sd(predictor, na.rm = TRUE)) || stats::sd(predictor, na.rm = TRUE) == 0) {
    return(list(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = length(response)))
  }

  df <- data.frame(y = response, predictor = predictor)
  if (!is.null(covariates) && ncol(covariates) > 0) {
    cov_df <- as.data.frame(covariates, stringsAsFactors = FALSE)
    colnames(cov_df) <- safe_make_names(colnames(cov_df))
    df <- cbind(df, cov_df)
  }

  fit <- tryCatch(stats::lm(y ~ ., data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(list(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(df)))
  }

  coef_tab <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
  if (is.null(coef_tab) || !("predictor" %in% rownames(coef_tab))) {
    return(list(beta = NA_real_, se = NA_real_, p_value = NA_real_, n = nrow(df)))
  }

  list(
    beta = unname(coef_tab["predictor", "Estimate"]),
    se = unname(coef_tab["predictor", "Std. Error"]),
    p_value = unname(coef_tab["predictor", "Pr(>|t|)"]),
    n = nrow(df),
    model = fit,
    term_name = predictor_name
  )
}

make_peak_table <- function(peak_coords, loci_df, frac_open_mat, mean_acc_mat) {
  peak_mid <- vapply(strsplit(peak_coords, "-"), function(parts) {
    as.integer((as.numeric(parts[2]) + as.numeric(parts[3])) %/% 2)
  }, FUN.VALUE = integer(1))
  peak_chr <- vapply(strsplit(peak_coords, "-"), function(parts) parts[1], FUN.VALUE = character(1))
  peak_start <- vapply(strsplit(peak_coords, "-"), function(parts) as.integer(parts[2]), FUN.VALUE = integer(1))
  peak_end <- vapply(strsplit(peak_coords, "-"), function(parts) as.integer(parts[3]), FUN.VALUE = integer(1))

  dist_to_locus <- sapply(seq_len(nrow(loci_df)), function(i) {
    ifelse(peak_chr == loci_df$chr[i], abs(peak_mid - loci_df$pos[i]), NA_integer_)
  })
  if (!is.matrix(dist_to_locus)) {
    dist_to_locus <- matrix(dist_to_locus, ncol = 1)
  }
  colnames(dist_to_locus) <- loci_df$rsid

  min_dist <- apply(dist_to_locus, 1, function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) NA_integer_ else min(v)
  })
  closest_idx <- apply(dist_to_locus, 1, function(v) {
    if (all(!is.finite(v))) return(NA_integer_)
    which.min(v)
  })
  closest_rsid <- ifelse(is.na(closest_idx), NA_character_, loci_df$rsid[closest_idx])

  candidate_table <- data.frame(
    peak = peak_coords,
    chr = peak_chr,
    start = peak_start,
    end = peak_end,
    mean_open_fraction = rowMeans(frac_open_mat, na.rm = TRUE),
    median_open_fraction = apply(frac_open_mat, 1, stats::median, na.rm = TRUE),
    mean_accessibility = rowMeans(mean_acc_mat, na.rm = TRUE),
    min_distance_to_locus = min_dist,
    closest_locus = closest_rsid,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(loci_df))) {
    candidate_table[[paste0("distance_", loci_df$rsid[i])]] <- dist_to_locus[, i]
  }
  candidate_table
}

run_phase_locus_tf_pipeline <- function(target_gene,
                                        loci_csv,
                                        phase,
                                        cell_type,
                                        window_bp = 5e4,
                                        enhancer_prox_bp = 5000,
                                        min_open_fraction = 0.01,
                                        fimo_hit_p_value = 1e-4,
                                        fimo_scan_p_value = 0.05,
                                        max_candidate_peaks = 50,
                                        n_cores = 10,
                                        n_runs = 20,
                                        p_value_cutoff = 0.05,
                                        output_root = RESULTS_ROOT) {
  addArchRThreads(threads = min(4L, max(1L, as.integer(n_cores))))
  addArchRGenome("hg38")

  phase <- tolower(phase)
  if (!(phase %in% c("prenatal", "postnatal"))) {
    stop("Phase must be prenatal or postnatal")
  }

  cell_type <- if (toupper(cell_type) == "ALL") "ALL" else c(
    "EXN" = "ExN",
    "INN" = "InN",
    "GLIA" = "Glia",
    "PROGENITOR" = "Progenitor",
    "NNC" = "NNC"
  )[toupper(cell_type)]
  if (is.na(cell_type)) {
    stop("CELL_TYPE must be one of ALL, ExN, InN, Glia, Progenitor, NNC")
  }

  loci_path <- if (grepl("^/", loci_csv)) loci_csv else file.path(PROJECT_ROOT, loci_csv)
  loci_df <- read_loci(loci_path)
  locus_tag <- tools::file_path_sans_ext(basename(loci_path))

  out_dir <- file.path(
    output_root,
    paste0("celltype_", cell_type),
    paste0("phase_", phase),
    target_gene,
    locus_tag
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  fimo_dir <- file.path(out_dir, "fimo_results")
  dir.create(fimo_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n", paste(rep("=", 88), collapse = ""), "\n", sep = "")
  cat("PHASE-LEVEL LOCUS TF PIPELINE\n")
  cat(paste(rep("=", 88), collapse = ""), "\n")
  cat("Target gene:", target_gene, "\n")
  cat("Phase:", phase, "\n")
  cat("Cell type:", cell_type, "\n")
  cat("Window around loci:", window_bp, "bp\n")
  cat("Enhancer proximity threshold:", enhancer_prox_bp, "bp\n")
  cat("Minimum open fraction:", min_open_fraction, "\n")
  cat("FIMO hit p-value threshold:", fimo_hit_p_value, "\n")
  cat("FIMO scan p-value threshold:", fimo_scan_p_value, "\n")
  cat("Max candidate peaks:", max_candidate_peaks, "\n")

  if (!file.exists(MAIN_ARCHR_PROJ_FILE)) {
    stop("Main ArchR project not found: ", MAIN_ARCHR_PROJ_FILE)
  }
  if (!dir.exists(JASPAR_DIR)) {
    stop("JASPAR directory not found: ", JASPAR_DIR)
  }

  batch_table <- get_phase_batch_table(phase)
  if (nrow(batch_table) == 0) {
    stop("No valid ATAC batches found for phase ", phase)
  }
  phase_periods <- get_phase_periods(phase)
  allowed_cell_ids <- get_archr_cells_for_phase(phase_periods, batch_table, cell_type)
  if (!is.null(allowed_cell_ids) && length(allowed_cell_ids) == 0) {
    cat("Warning: no ATAC metadata keys found for selected cell type; falling back to mapping-only matching.\n")
    allowed_cell_ids <- NULL
  }

  archr_proj <- readRDS(MAIN_ARCHR_PROJ_FILE)
  peak_set <- getPeakSet(archr_proj)
  if (is.null(peak_set) || length(peak_set) == 0) {
    stop("Peak set missing from main ArchR project")
  }

  loci_ranges <- GRanges(
    seqnames = loci_df$chr,
    ranges = IRanges(start = pmax(1L, loci_df$pos - as.integer(window_bp)),
                     end = loci_df$pos + as.integer(window_bp)),
    rsid = loci_df$rsid
  )
  ov <- findOverlaps(peak_set, loci_ranges)
  peak_idx <- unique(queryHits(ov))
  if (length(peak_idx) == 0) {
    stop("No ArchR peaks overlap the supplied loci windows")
  }

  candidate_peak_ranges <- peak_set[peak_idx]
  candidate_peak_coords <- paste(seqnames(candidate_peak_ranges), start(candidate_peak_ranges), end(candidate_peak_ranges), sep = "-")
  candidate_peak_norm <- normalize_peak_id(candidate_peak_coords)
  batch_keys <- paste0("p", batch_table$period, "::", batch_table$batch)

  frac_open_mat <- matrix(
    NA_real_,
    nrow = length(candidate_peak_coords),
    ncol = nrow(batch_table),
    dimnames = list(candidate_peak_coords, batch_keys)
  )
  mean_acc_mat <- matrix(
    NA_real_,
    nrow = length(candidate_peak_coords),
    ncol = nrow(batch_table),
    dimnames = list(candidate_peak_coords, batch_keys)
  )

  cat("Computing phase-level peak openness across", nrow(batch_table), "period+batch groups...\n")
  for (i in seq_len(nrow(batch_table))) {
    p <- batch_table$period[i]
    b <- batch_table$batch[i]
    key <- batch_keys[i]

    proj_batch <- tryCatch(helpers$load_batch_archr_from_arrow(b), error = function(e) NULL)
    if (is.null(proj_batch)) {
      next
    }
    peak_se <- tryCatch(getMatrixFromProject(ArchRProj = proj_batch, useMatrix = "PeakMatrix"), error = function(e) NULL)
    if (is.null(peak_se)) {
      next
    }
    peak_mat <- assay(peak_se)
    if (ncol(peak_mat) == 0) {
      next
    }

    row_norm <- get_peak_ids_from_se(peak_se)
    row_idx <- match(candidate_peak_norm, row_norm)
    valid_rows <- !is.na(row_idx)
    if (!any(valid_rows)) {
      next
    }

    batch_cells <- colnames(peak_mat)
    if (!is.null(allowed_cell_ids)) {
      batch_cell_ids <- helpers$canonical_cell_id(batch_cells, sample_hint = b)
      keep_allowed <- batch_cell_ids %in% allowed_cell_ids
      if (any(keep_allowed)) {
        batch_cells <- batch_cells[keep_allowed]
      }
    }
    if (length(batch_cells) == 0) {
      next
    }

    mc_groups <- build_batch_metacell_groups(batch_name = b, atac_cells = batch_cells, cell_type = cell_type)
    mc_groups <- mc_groups[lengths(mc_groups) > 0]
    if (length(mc_groups) == 0) {
      next
    }

    sub_mat <- peak_mat[row_idx[valid_rows], , drop = FALSE]
    open_bin <- matrix(0, nrow = sum(valid_rows), ncol = length(mc_groups))
    mean_mc <- matrix(NA_real_, nrow = sum(valid_rows), ncol = length(mc_groups))

    j <- 0L
    for (mc_cells in mc_groups) {
      idx <- match(mc_cells, colnames(peak_mat))
      idx <- idx[!is.na(idx)]
      if (length(idx) == 0) {
        next
      }
      j <- j + 1L
      cell_mat <- sub_mat[, idx, drop = FALSE]
      open_bin[, j] <- as.numeric(Matrix::rowSums(cell_mat > 0) > 0)
      mean_mc[, j] <- as.numeric(Matrix::rowMeans(cell_mat))
    }
    if (j == 0L) {
      next
    }

    frac_open_mat[valid_rows, key] <- rowMeans(open_bin[, seq_len(j), drop = FALSE], na.rm = TRUE)
    mean_acc_mat[valid_rows, key] <- rowMeans(mean_mc[, seq_len(j), drop = FALSE], na.rm = TRUE)
  }

  if (all(!is.finite(frac_open_mat))) {
    stop("No valid batch/cell values to compute phase-specific peak openness")
  }

  candidate_peak_table <- make_peak_table(candidate_peak_coords, loci_df, frac_open_mat, mean_acc_mat)
  candidate_peak_table <- candidate_peak_table[order(candidate_peak_table$min_distance_to_locus, -candidate_peak_table$mean_open_fraction), , drop = FALSE]
  write.csv(candidate_peak_table, file.path(out_dir, paste0(target_gene, "_candidate_peaks.csv")), row.names = FALSE)

  priority_keep <- with(
    candidate_peak_table,
    is.finite(min_distance_to_locus) &
      min_distance_to_locus <= enhancer_prox_bp &
      is.finite(mean_open_fraction) &
      mean_open_fraction >= min_open_fraction
  )
  if (!any(priority_keep)) {
    fallback_idx <- which(is.finite(candidate_peak_table$min_distance_to_locus))
    fallback_idx <- fallback_idx[order(candidate_peak_table$min_distance_to_locus[fallback_idx], -candidate_peak_table$mean_open_fraction[fallback_idx])]
    fallback_idx <- head(fallback_idx, min(20L, length(fallback_idx)))
    priority_peak_table <- candidate_peak_table[fallback_idx, , drop = FALSE]
  } else {
    priority_peak_table <- candidate_peak_table[priority_keep, , drop = FALSE]
    priority_peak_table <- priority_peak_table[order(priority_peak_table$min_distance_to_locus, -priority_peak_table$mean_open_fraction), , drop = FALSE]
  }

  if (nrow(priority_peak_table) == 0) {
    stop("No priority peaks remained after locus/accessibility filtering")
  }
  if (nrow(priority_peak_table) > max_candidate_peaks) {
    priority_peak_table <- head(priority_peak_table, max_candidate_peaks)
  }
  write.csv(priority_peak_table, file.path(out_dir, paste0(target_gene, "_priority_peaks.csv")), row.names = FALSE)

  cat("Candidate peaks in loci windows:", nrow(candidate_peak_table), "\n")
  cat("Priority peaks for motif scan/regression:", nrow(priority_peak_table), "\n")

  peak_gr <- GRanges(
    seqnames = priority_peak_table$chr,
    ranges = IRanges(priority_peak_table$start, priority_peak_table$end)
  )
  peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, peak_gr)
  names(peak_seqs) <- priority_peak_table$peak
  fasta_file <- file.path(out_dir, paste0(target_gene, "_priority_peaks.fa"))
  writeXStringSet(peak_seqs, fasta_file)

  bg_file <- file.path(out_dir, paste0(target_gene, "_priority_peaks.bg"))
  bg_built <- helpers$build_markov_bgfile(fasta_filepath = fasta_file, bgfile_path = bg_file, order = 1L)

  tf_collection <- helpers$extract_all_tf_names(JASPAR_DIR)
  motif_library <- vector("list", length(tf_collection))
  names(motif_library) <- tf_collection

  cat("Running FIMO for", length(tf_collection), "TF motifs...\n")
  run_one_tf <- function(tf_name) {
    bvec <- tryCatch(
      helpers$run_fimo_pipeline_optimized(
        tf_folder = JASPAR_DIR,
        tf_name = tf_name,
        fimo_folder = fimo_dir,
        coords = priority_peak_table$peak,
        tg_name = target_gene,
        p_value_threshold = fimo_hit_p_value,
        fasta_filepath = fasta_file,
        scan_p_value_threshold = fimo_scan_p_value,
        bgfile_path = bg_built
      ),
      error = function(e) NULL
    )
    if (is.null(bvec)) {
      bvec <- rep(FALSE, nrow(priority_peak_table))
    }
    as.integer(bvec)
  }

  max_cores <- max(1L, min(as.integer(n_cores), parallel::detectCores(logical = TRUE), length(tf_collection)))
  if (max_cores > 1L && .Platform$OS.type != "windows") {
    motif_res <- parallel::mclapply(tf_collection, run_one_tf, mc.cores = max_cores, mc.preschedule = FALSE)
    motif_library <- motif_res
  } else {
    for (i in seq_along(tf_collection)) {
      motif_library[[i]] <- run_one_tf(tf_collection[i])
      if (i %% 100 == 0 || i == length(tf_collection)) {
        cat("  FIMO progress:", i, "/", length(tf_collection), "\n")
      }
    }
  }

  motif_matrix <- do.call(rbind, motif_library)
  rownames(motif_matrix) <- tf_collection
  colnames(motif_matrix) <- priority_peak_table$peak
  storage.mode(motif_matrix) <- "integer"
  motif_matrix <- motif_matrix[rowSums(motif_matrix) > 0, , drop = FALSE]
  if (nrow(motif_matrix) == 0) {
    stop("No TF motifs matched the selected priority peaks")
  }
  write.csv(motif_matrix, file.path(out_dir, paste0(target_gene, "_motif_matrix.csv")))

  tf_names <- rownames(motif_matrix)
  gene_request <- unique(c(tf_names, target_gene))

  X_list <- list()
  Y_list <- list()
  Z_list <- list()
  batch_info_rows <- list()

  cat("Retrieving metacell RNA and ATAC matrices across phase...\n")
  for (p in phase_periods) {
    batches <- get_period_batches(p)
    if (length(batches) == 0) {
      next
    }

    expr_all <- helpers$load_rna_expression_for_period(p, gene_request, temp_dir = TMP_ROOT)
    if (is.null(expr_all)) {
      cat("  Warning: RNA extraction failed for period", p, "\n")
      next
    }
    missing_genes <- setdiff(gene_request, colnames(expr_all))
    if (length(missing_genes) > 0) {
      present_tfs <- intersect(tf_names, colnames(expr_all))
    } else {
      present_tfs <- tf_names
    }
    if (!(target_gene %in% colnames(expr_all))) {
      cat("  Warning: target gene missing from RNA for period", p, "\n")
      next
    }
    if (length(present_tfs) == 0) {
      cat("  Warning: no motif-hit TFs present in RNA for period", p, "\n")
      next
    }

    X_cells <- expr_all[, present_tfs, drop = FALSE]
    Y_cells <- expr_all[, target_gene, drop = FALSE]

    for (b in batches) {
      mc_assign <- helpers$load_metacell_assignments_filtered(b, cell_type = cell_type)
      if (is.null(mc_assign) || length(mc_assign) == 0) {
        next
      }

      X_batch <- helpers$aggregate_rna_to_metacells(X_cells, mc_assign)
      Y_batch <- helpers$aggregate_rna_to_metacells(Y_cells, mc_assign)
      proj_batch <- tryCatch(helpers$load_batch_archr_from_arrow(b), error = function(e) NULL)
      if (is.null(proj_batch)) {
        next
      }
      Z_batch <- tryCatch(
        helpers$aggregate_atac_to_metacells(proj_batch, b, priority_peak_table$peak, mc_assign),
        error = function(e) NULL
      )
      if (is.null(Z_batch)) {
        next
      }

      row_id <- paste0("p", p, "::", b, "::", rownames(X_batch))
      rownames(X_batch) <- row_id
      rownames(Y_batch) <- row_id
      rownames(Z_batch) <- row_id

      X_list[[paste0("p", p, "::", b)]] <- X_batch
      Y_list[[paste0("p", p, "::", b)]] <- Y_batch
      Z_list[[paste0("p", p, "::", b)]] <- Z_batch
      batch_info_rows[[paste0("p", p, "::", b)]] <- data.frame(
        metacell_id = row_id,
        period = as.integer(p),
        batch = b,
        batch_phase = paste0("p", p, "::", b),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(X_list) == 0 || length(Y_list) == 0 || length(Z_list) == 0) {
    stop("No valid phase-level metacell matrices could be built for this target/phase/cell type")
  }

  common_tfs <- Reduce(intersect, lapply(X_list, colnames))
  common_peaks <- Reduce(intersect, lapply(Z_list, colnames))
  if (length(common_tfs) == 0) {
    stop("No common TF columns across phase batches")
  }
  if (length(common_peaks) == 0) {
    stop("No common peak columns across phase batches")
  }

  X_raw <- do.call(rbind, lapply(X_list, function(x) x[, common_tfs, drop = FALSE]))
  Y_raw <- do.call(rbind, Y_list)
  Z_raw <- do.call(rbind, lapply(Z_list, function(z) z[, common_peaks, drop = FALSE]))
  batch_info <- do.call(rbind, batch_info_rows)
  batch_info <- batch_info[match(rownames(X_raw), batch_info$metacell_id), , drop = FALSE]

  x_upper <- toupper(colnames(X_raw))
  target_upper <- toupper(target_gene)
  target_idx <- which(x_upper == target_upper)
  if (length(target_idx) > 0) {
    X_raw <- X_raw[, -target_idx, drop = FALSE]
  }

  n_batch_groups <- length(unique(batch_info$batch_phase))
  if (n_batch_groups > 1) {
    tmp_target_col <- "TARGETGENEFORRPCA"
    Y_for_correction <- Y_raw
    colnames(Y_for_correction) <- tmp_target_col
    XY_raw <- cbind(X_raw, Y_for_correction)
    XY_corrected <- tryCatch(
      helpers$apply_rna_batch_correction(XY_raw, batch_info$batch_phase),
      error = function(e) {
        cat("Warning: RNA batch correction failed; using uncorrected RNA. Reason:", conditionMessage(e), "\n")
        XY_raw
      }
    )
    X_corrected <- XY_corrected[, colnames(X_raw), drop = FALSE]
    Y_corrected <- XY_corrected[, tmp_target_col, drop = FALSE]
    colnames(Y_corrected) <- target_gene
    Z_corrected <- tryCatch(
      helpers$apply_atac_batch_correction(as.matrix(Z_raw), batch_info$batch_phase),
      error = function(e) {
        cat("Warning: ATAC batch correction failed; using uncorrected ATAC. Reason:", conditionMessage(e), "\n")
        as.matrix(Z_raw)
      }
    )
  } else {
    X_corrected <- as.matrix(X_raw)
    Y_corrected <- as.matrix(Y_raw)
    Z_corrected <- as.matrix(Z_raw)
  }

  X_scaled <- scale_matrix_columns(X_corrected, log1p_nonnegative = FALSE)
  Y_scaled <- scale_matrix_columns(Y_corrected, log1p_nonnegative = FALSE)
  Z_scaled <- scale_matrix_columns(Z_corrected, log1p_nonnegative = TRUE)

  write.csv(X_corrected, file.path(out_dir, paste0(target_gene, "_X_matrix_raw.csv")))
  write.csv(Y_corrected, file.path(out_dir, paste0(target_gene, "_Y_matrix_raw.csv")))
  write.csv(Z_corrected, file.path(out_dir, paste0(target_gene, "_Z_matrix_raw.csv")))
  write.csv(X_scaled, file.path(out_dir, paste0(target_gene, "_X_matrix_scaled.csv")))
  write.csv(Y_scaled, file.path(out_dir, paste0(target_gene, "_Y_matrix_scaled.csv")))
  write.csv(Z_scaled, file.path(out_dir, paste0(target_gene, "_Z_matrix_scaled.csv")))
  write.csv(batch_info, file.path(out_dir, paste0(target_gene, "_batch_info.csv")), row.names = FALSE)

  aligned <- helpers$verify_matrices(list(
    X_matrix = X_scaled,
    Y_matrix = Y_scaled,
    Z_matrix = Z_scaled,
    motif_matrix = motif_matrix,
    batch_info = batch_info
  ))

  period_covariates <- NULL
  if (!is.null(aligned$batch_info) && "period" %in% colnames(aligned$batch_info)) {
    period_covariates <- stats::model.matrix(~ 0 + factor(period), data = aligned$batch_info)
    colnames(period_covariates) <- paste0("period_", sub("^factor\\(period\\)", "", colnames(period_covariates)))
  }

  progress_file <- file.path(out_dir, paste0(target_gene, "_bootstrap_progress.txt"))
  cat("Running GRN-style inference on locus-restricted matrices...\n")
  grn_results <- run_grn_inference_with_scores_locus(
    helpers = helpers,
    X_matrix = aligned$X_matrix,
    Y_matrix = aligned$Y_matrix,
    Z_matrix = aligned$Z_matrix,
    motif_matrix = aligned$motif_matrix,
    batch_covariates = period_covariates,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff,
    progress_file = progress_file
  )

  if (!is.null(grn_results$bootstrap_details)) {
    write.csv(grn_results$bootstrap_details,
              file.path(out_dir, paste0(target_gene, "_bootstrap_variable_details.csv")),
              row.names = FALSE)
  }

  if (!is.null(grn_results$score_summary)) {
    write.csv(grn_results$score_summary,
              file.path(out_dir, paste0(target_gene, "_variable_score_summary.csv")),
              row.names = FALSE)
  }

  sig_tfs_only_file <- file.path(out_dir, paste0(target_gene, "_significant_tfs.csv"))
  if (nrow(grn_results$significant_tfs_only_df) > 0) {
    sig_tfs_only_df <- grn_results$significant_tfs_only_df[, c("Variable", "Mean_Coefficient", "Median_Coefficient", "Selection_Frequency", "Adjusted_P_value"), drop = FALSE]
    colnames(sig_tfs_only_df) <- c("TF", "Mean_Score", "Median_Score", "Selection_Frequency", "P_value")
    write.csv(sig_tfs_only_df, sig_tfs_only_file, row.names = FALSE)
  } else {
    write.csv(data.frame(TF = character(0), P_value = numeric(0), stringsAsFactors = FALSE),
              sig_tfs_only_file, row.names = FALSE)
  }

  sig_cres_file <- file.path(out_dir, paste0(target_gene, "_significant_cres.csv"))
  if (nrow(grn_results$significant_cres_df) > 0) {
    sig_cres_df <- grn_results$significant_cres_df[, c("Variable", "Mean_Coefficient", "Median_Coefficient", "Selection_Frequency", "Adjusted_P_value"), drop = FALSE]
    colnames(sig_cres_df) <- c("CRE", "Mean_Score", "Median_Score", "Selection_Frequency", "P_value")
    peak_map <- priority_peak_table[match(sig_cres_df$CRE, priority_peak_table$peak), , drop = FALSE]
    sig_cres_df$Closest_Locus <- peak_map$closest_locus
    sig_cres_df$Min_Distance_To_Locus <- peak_map$min_distance_to_locus
    sig_cres_df$Mean_Open_Fraction <- peak_map$mean_open_fraction
    write.csv(sig_cres_df, sig_cres_file, row.names = FALSE)
  } else {
    write.csv(data.frame(CRE = character(0), P_value = numeric(0), stringsAsFactors = FALSE),
              sig_cres_file, row.names = FALSE)
  }

  identified_tfs_file <- file.path(out_dir, paste0(target_gene, "_identified_tfs.csv"))
  if (nrow(grn_results$identified_tfs_df) > 0) {
    identified_tfs_df <- grn_results$identified_tfs_df[, c("Variable", "Mean_Coefficient", "Median_Coefficient", "Selection_Frequency", "Adjusted_P_value"), drop = FALSE]
    colnames(identified_tfs_df) <- c("TF", "Mean_Score", "Median_Score", "Selection_Frequency", "P_value")
    write.csv(identified_tfs_df, identified_tfs_file, row.names = FALSE)
  } else {
    write.csv(data.frame(TF = character(0), P_value = numeric(0), stringsAsFactors = FALSE),
              identified_tfs_file, row.names = FALSE)
  }

  combined_pval_df <- data.frame(
    Variable = names(grn_results$combined_p_values),
    P_value = grn_results$combined_p_values,
    Adjusted_P_value = grn_results$adjusted_p_values,
    stringsAsFactors = FALSE
  )
  write.csv(combined_pval_df, file.path(out_dir, paste0(target_gene, "_combined_p_values.csv")), row.names = FALSE)

  if (!is.null(grn_results$bootstrap_diagnostics)) {
    write.csv(grn_results$bootstrap_diagnostics,
              file.path(out_dir, paste0(target_gene, "_bootstrap_diagnostics.csv")),
              row.names = FALSE)
  }

  metadata <- list(
    target_gene = target_gene,
    phase = phase,
    cell_type = cell_type,
    loci_csv = loci_path,
    loci = loci_df,
    window_bp = window_bp,
    enhancer_prox_bp = enhancer_prox_bp,
    min_open_fraction = min_open_fraction,
    fimo_hit_p_value = fimo_hit_p_value,
    fimo_scan_p_value = fimo_scan_p_value,
    max_candidate_peaks = max_candidate_peaks,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff,
    n_candidate_peaks = nrow(candidate_peak_table),
    n_priority_peaks = nrow(priority_peak_table),
    n_tfs_with_motif = nrow(motif_matrix),
    n_metacells = nrow(aligned$X_matrix),
    n_significant_tfs = ifelse(is.null(grn_results$n_sig_tfs), length(grn_results$significant_tfs_only), grn_results$n_sig_tfs),
    n_significant_cres = ifelse(is.null(grn_results$n_sig_cres), length(grn_results$significant_cres), grn_results$n_sig_cres),
    n_identified_tfs = ifelse(is.null(grn_results$n_identified_tfs), length(grn_results$significant_tfs), grn_results$n_identified_tfs)
  )
  saveRDS(metadata, file.path(out_dir, paste0(target_gene, "_metadata.rds")))

  summary_lines <- c(
    "Phase-Level Locus TF Pipeline Summary",
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Phase:", phase),
    paste("Cell Type:", cell_type),
    paste("Loci:", paste(loci_df$rsid, collapse = ", ")),
    paste("Locus Window:", window_bp, "bp"),
    paste("Enhancer Proximity Threshold:", enhancer_prox_bp, "bp"),
    paste("Peak Open Fraction Threshold:", min_open_fraction),
    paste("FIMO Hit P-value Threshold:", fimo_hit_p_value),
    paste("FIMO Scan P-value Threshold:", fimo_scan_p_value),
    paste("Bootstrap Runs:", n_runs),
    paste("P-value Cutoff (BH):", p_value_cutoff),
    "",
    "Normalization and regression:",
    "  X: phase-level RNA batch correction (Seurat RPCA when >1 batch), then column z-score",
    "  Z: phase-level ATAC batch correction (LSI when >1 batch), log1p transform, then column z-score",
    "  Y: column z-score",
    "  Regression: overlapped group lasso + post-selection inference",
    "  Period covariates: included as unpenalized variables",
    "  P-values: combined across bootstrap runs with Cauchy combination",
    "  Scores: mean/median coefficient across selected bootstrap runs only",
    "",
    paste("Candidate peaks:", nrow(candidate_peak_table)),
    paste("Priority peaks:", nrow(priority_peak_table)),
    paste("Motif-supported TFs:", nrow(motif_matrix)),
    paste("Metacells used:", nrow(aligned$X_matrix)),
    paste("Significant TFs (p-value < cutoff):", metadata$n_significant_tfs),
    paste("Significant CREs (p-value < cutoff):", metadata$n_significant_cres),
    paste("Identified TFs (TF and binding CRE significant):", metadata$n_identified_tfs),
    "",
    paste("Output Directory:", out_dir)
  )
  writeLines(summary_lines, file.path(out_dir, paste0(target_gene, "_summary.txt")))

  cat("\nPipeline complete.\n")
  cat("Output directory:", out_dir, "\n")

  invisible(list(
    candidate_peaks = candidate_peak_table,
    priority_peaks = priority_peak_table,
    motif_matrix = motif_matrix,
    X_scaled = X_scaled,
    Y_scaled = Y_scaled,
    Z_scaled = Z_scaled,
    grn_results = grn_results,
    metadata = metadata
  ))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  if (length(args) < 4) {
    stop(
      "Usage: Rscript coloc/run_locus_tf.R TARGET_GENE LOCI_CSV PHASE CELL_TYPE ",
      "[WINDOW_BP] [ENHANCER_PROX_BP] [MIN_OPEN_FRACTION] [FIMO_HIT_P] [FIMO_SCAN_P] ",
      "[MAX_CANDIDATE_PEAKS] [N_CORES] [N_RUNS] [P_VALUE_CUTOFF]"
    )
  }

  target_gene <- args[1]
  loci_csv <- args[2]
  phase <- args[3]
  cell_type <- args[4]
  window_bp <- if (length(args) >= 5) as.numeric(args[5]) else 5e4
  enhancer_prox_bp <- if (length(args) >= 6) as.numeric(args[6]) else 5000
  min_open_fraction <- if (length(args) >= 7) as.numeric(args[7]) else 0.01
  fimo_hit_p_value <- if (length(args) >= 8) as.numeric(args[8]) else 1e-4
  fimo_scan_p_value <- if (length(args) >= 9) as.numeric(args[9]) else 0.05
  max_candidate_peaks <- if (length(args) >= 10) as.integer(args[10]) else 50L
  n_cores <- if (length(args) >= 11) as.integer(args[11]) else 10L
  n_runs <- if (length(args) >= 12) as.integer(args[12]) else 20L
  p_value_cutoff <- if (length(args) >= 13) as.numeric(args[13]) else 0.05

  run_phase_locus_tf_pipeline(
    target_gene = target_gene,
    loci_csv = loci_csv,
    phase = phase,
    cell_type = cell_type,
    window_bp = window_bp,
    enhancer_prox_bp = enhancer_prox_bp,
    min_open_fraction = min_open_fraction,
    fimo_hit_p_value = fimo_hit_p_value,
    fimo_scan_p_value = fimo_scan_p_value,
    max_candidate_peaks = max_candidate_peaks,
    n_cores = n_cores,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff
  )
} else {
  cat(
    "Usage: Rscript coloc/run_locus_tf.R TARGET_GENE LOCI_CSV PHASE CELL_TYPE ",
    "[WINDOW_BP] [ENHANCER_PROX_BP] [MIN_OPEN_FRACTION] [FIMO_HIT_P] [FIMO_SCAN_P] ",
    "[MAX_CANDIDATE_PEAKS] [N_CORES] [N_RUNS] [P_VALUE_CUTOFF]\n",
    sep = ""
  )
}
