#!/usr/bin/env Rscript

# =============================================================================
# DATA RETRIEVAL HELPERS
# =============================================================================
# This script retrieves RNA expression (X, Y) and ATAC accessibility (Z) data
# at the metacell level for all batches in a period, then performs batch correction
# and aggregates the results.
#
# Internal module shared by the phase and locus entry points.

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(GenomicRanges)
  library(Matrix)
  library(Seurat)
  library(irlba)
})

# Initialize ArchR only if not being sourced
if (!exists(".SOURCING_PIPELINES")) {
  addArchRThreads(threads = 4)
  addArchRGenome("hg38")
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths are supplied by the public entry point when this module is sourced.
if (!exists("BASE_DIR", inherits = FALSE)) BASE_DIR <- getwd()
if (!exists("DATA_DIR", inherits = FALSE)) DATA_DIR <- file.path(BASE_DIR, "../data")
if (!exists("ATAC_DIR", inherits = FALSE)) ATAC_DIR <- file.path(DATA_DIR, "atac")
if (!exists("PERIOD_MOTIF_RESULTS_DIR", inherits = FALSE)) PERIOD_MOTIF_RESULTS_DIR <- file.path(BASE_DIR, "period_motif_matching_results")
if (!exists("METACELL_DIR", inherits = FALSE)) METACELL_DIR <- file.path(BASE_DIR, "../explore/explore_results/metacells")
if (!exists("OUTPUT_DIR", inherits = FALSE)) OUTPUT_DIR <- file.path(BASE_DIR, "period_data_retrieval_results")

# Main ArchRProject file
if (!exists("MAIN_ARCHR_PROJ_FILE", inherits = FALSE)) MAIN_ARCHR_PROJ_FILE <- file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-Save-ArchR-Project.rds")
if (!exists("ATAC_META_FILE", inherits = FALSE)) ATAC_META_FILE <- file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-ATAC_P3-P15_meta_D240727.tsv")

# Python helper script
PYTHON_EXTRACT_RNA_SCRIPT <- file.path(
  if (exists("INTERNAL_R_DIR", inherits = FALSE)) INTERNAL_R_DIR else BASE_DIR,
  "extract_rna_expression.py"
)

.MAIN_PEAK_SET_CACHE <- NULL

get_main_peak_set_cached <- function() {
  if (!is.null(.MAIN_PEAK_SET_CACHE)) {
    return(.MAIN_PEAK_SET_CACHE)
  }
  if (!file.exists(MAIN_ARCHR_PROJ_FILE)) {
    return(NULL)
  }
  main_proj <- tryCatch(readRDS(MAIN_ARCHR_PROJ_FILE), error = function(e) NULL)
  if (is.null(main_proj)) return(NULL)
  peak_set <- tryCatch(getPeakSet(main_proj), error = function(e) NULL)
  if (!is.null(peak_set) && length(peak_set) > 0) {
    .MAIN_PEAK_SET_CACHE <<- peak_set
  }
  .MAIN_PEAK_SET_CACHE
}

compact_signature <- function(x) {
  x <- as.character(x)
  n <- length(x)
  if (n == 0) return("n=0")
  idx <- unique(round(seq(1, n, length.out = min(20, n))))
  paste0("n=", n, ";sample=", paste(x[idx], collapse = "|"))
}

# =============================================================================
# HELPER FUNCTIONS (from data_retrieval_pipeline.R)
# =============================================================================

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

extract_barcode <- function(x) {
  x <- as.character(x)
  x <- sub("^.*[_#]", "", x)
  sub("-1$", "", x)
}

extract_sample <- function(x) {
  x <- as.character(x)
  out <- ifelse(grepl("#", x), sub("#.*$", "", x), ifelse(grepl("_", x), sub("_.*$", "", x), NA_character_))
  out[!nzchar(out)] <- NA_character_
  out
}

canonical_cell_id <- function(x, sample_hint = NA_character_) {
  x <- as.character(x)
  sample <- extract_sample(x)
  barcode <- extract_barcode(x)
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

#' Get all batches belonging to a specific period
#' @param period Period number (e.g., 3, 4, 5)
#' @return Vector of batch names
get_valid_atac_batches_for_period <- function(period) {
  if (!file.exists(ATAC_META_FILE)) {
    stop("ERROR: ATAC metadata file not found: ", ATAC_META_FILE)
  }

  atac_meta <- tryCatch({
    read_atac_meta(ATAC_META_FILE)
  }, error = function(e) {
    stop("ERROR: Failed to read ATAC metadata file ", ATAC_META_FILE, ": ", e$message)
  })

  if (!("Sample" %in% colnames(atac_meta)) || !("Period" %in% colnames(atac_meta))) {
    stop("ERROR: ATAC metadata must contain 'Sample' and 'Period' columns: ", ATAC_META_FILE)
  }

  period_num <- as.numeric(period)
  meta_period <- suppressWarnings(as.numeric(as.character(atac_meta$Period)))
  period_batches <- unique(as.character(atac_meta$Sample[!is.na(meta_period) & meta_period == period_num]))
  period_batches <- period_batches[nzchar(period_batches)]

  arrow_pattern <- "ArrowFiles-.*\\.arrow$"
  arrow_files <- list.files(ATAC_DIR, pattern = arrow_pattern, full.names = FALSE)
  arrow_batches <- sub("\\.arrow$", "", sub("^.*ArrowFiles-", "", arrow_files))

  sort(intersect(period_batches, arrow_batches))
}

get_batches_for_period <- function(period) {
  pattern <- paste0("^mapping_p", period, "_.*\\.csv$")
  mapping_files <- list.files(METACELL_DIR, pattern = pattern, full.names = FALSE)

  if (length(mapping_files) == 0) {
    return(character(0))
  }

  batch_names <- gsub(paste0("^mapping_p", period, "_"), "", mapping_files)
  batch_names <- gsub("\\.csv$", "", batch_names)

  valid_atac_batches <- get_valid_atac_batches_for_period(period)
  filtered_batches <- sort(intersect(batch_names, valid_atac_batches))
  removed_batches <- setdiff(batch_names, filtered_batches)

  if (length(removed_batches) > 0) {
    cat("Warning: Excluding mismatched/missing-Arrow batches for period", period, ":",
        paste(sort(removed_batches), collapse = ", "), "\n")
  }

  return(filtered_batches)
}

#' Check if period motif matching results exist
#' @param target_gene Target gene symbol
#' @param period Period number
#' @return List with exists status and file paths
check_period_motif_results <- function(target_gene, period, motif_results_dir = PERIOD_MOTIF_RESULTS_DIR) {
  period_dir <- file.path(motif_results_dir, paste0("period_", period), target_gene)
  motif_matrix_file <- file.path(period_dir, paste0(target_gene, "_motif_matrix.csv"))

  # Check if base directory exists
  if (!dir.exists(motif_results_dir)) {
    return(list(
      exists = FALSE,
      message = paste("Period motif matching results base directory not found:", motif_results_dir,
                     "\nCurrent working directory:", getwd(),
                     "\nBASE_DIR:", BASE_DIR)
    ))
  }

  # Check if period directory exists
  if (!dir.exists(period_dir)) {
    # List what's actually in the base directory to help debug
    period_subdirs <- list.dirs(motif_results_dir, full.names = FALSE, recursive = FALSE)
    target_subdirs <- character(0)
    if (length(period_subdirs) > 0) {
      # Check if period_3 directory exists
      period_pattern <- paste0("^period_", period, "$")
      matching_periods <- grep(period_pattern, period_subdirs, value = TRUE)
      if (length(matching_periods) > 0) {
        period_subdir_path <- file.path(motif_results_dir, matching_periods[1])
        target_subdirs <- list.dirs(period_subdir_path, full.names = FALSE, recursive = FALSE)
      }
    }

    return(list(
      exists = FALSE,
      message = paste("Period motif matching results directory not found:", period_dir,
                     "\nPERIOD_MOTIF_RESULTS_DIR:", motif_results_dir,
                     "\nAvailable period directories:", paste(period_subdirs, collapse = ", "),
                     ifelse(length(target_subdirs) > 0,
                            paste("\nTarget gene directories in period_", period, ":", paste(target_subdirs, collapse = ", ")),
                            ""))
    ))
  }

  # Check if motif matrix file exists
  if (!file.exists(motif_matrix_file)) {
    # List files in the directory to help debug
    files_in_dir <- list.files(period_dir, pattern = "\\.csv$")
    return(list(
      exists = FALSE,
      message = paste("Motif matrix file not found:", motif_matrix_file,
                     "\nFiles in directory:", paste(files_in_dir, collapse = ", "))
    ))
  }

  return(list(
    exists = TRUE,
    period_dir = period_dir,
    motif_matrix_file = motif_matrix_file,
    message = paste("Period motif matching results found for", target_gene, "in period", period)
  ))
}

#' Load RNA expression from h5ad file using Python helper
#' @param period Period number
#' @param gene_names Vector of gene names to extract
#' @param temp_dir Temporary directory for intermediate files
#' @return Expression matrix (metacells × genes) or NULL if failed
load_rna_expression_for_period <- function(period, gene_names, temp_dir = tempdir()) {
  h5ad_file <- file.path(METACELL_DIR, paste0("metacells_period_", period, ".h5ad"))

  if (!file.exists(h5ad_file)) {
    cat("    Error: h5ad file not found:", h5ad_file, "\n")
    cat("    METACELL_DIR:", METACELL_DIR, "\n")
    return(NULL)
  }

  if (!file.exists(PYTHON_EXTRACT_RNA_SCRIPT)) {
    cat("    Error: Python helper script not found:", PYTHON_EXTRACT_RNA_SCRIPT, "\n")
    return(NULL)
  }

  # Create temporary files with unique names to avoid conflicts when called multiple times
  temp_suffix <- paste0("_", Sys.getpid(), "_", sample(10000:99999, 1))
  gene_names_file <- file.path(temp_dir, paste0("genes_", period, temp_suffix, ".txt"))
  output_file <- file.path(temp_dir, paste0("expression_", period, temp_suffix, ".csv"))

  # Write gene names to file
  writeLines(gene_names, gene_names_file)

  # Run Python script and capture both stdout and stderr
  python_cmd <- sprintf("python %s %s %s %s 2>&1",
                        shQuote(PYTHON_EXTRACT_RNA_SCRIPT),
                        shQuote(h5ad_file),
                        shQuote(gene_names_file),
                        shQuote(output_file))

  result <- system(python_cmd, intern = TRUE)
  exit_status <- attr(result, "status")

  # Check for "No genes matched" error in output - this means gene not found
  if (!is.null(result) && length(result) > 0) {
    no_genes_matched <- any(grepl("ERROR.*No genes matched|No genes matched", result, ignore.case = TRUE))
    if (no_genes_matched) {
      # Extract the gene name(s) that were requested
      requested_genes <- if (length(gene_names) == 1) gene_names else paste(length(gene_names), "genes")
      error_msg <- paste0(
        "ERROR: Gene(s) not found in RNA expression data.\n",
        "    Requested: ", requested_genes, "\n",
        "    Period: ", period, "\n",
        "    The pipeline will stop to prevent hanging."
      )
      cat("    ", error_msg, "\n")
      # Clean up temp files
      unlink(c(gene_names_file, output_file))
      stop(error_msg)
    }

    # Print output for debugging
    cat("    Python script output:\n")
    for (line in result) {
      cat("    ", line, "\n")
    }
  }

  if (!is.null(exit_status) && exit_status != 0) {
    # Extract gene name(s) for error message
    requested_genes <- if (length(gene_names) == 1) gene_names else paste(length(gene_names), "genes")
    error_msg <- paste0(
      "ERROR: Python script failed to extract RNA expression.\n",
      "    Requested: ", requested_genes, "\n",
      "    Period: ", period, "\n",
      "    Exit status: ", exit_status, "\n",
      "    The pipeline will stop to prevent hanging."
    )
    cat("    ", error_msg, "\n")
    cat("    Command was:", python_cmd, "\n")
    # Clean up temp files
    unlink(c(gene_names_file, output_file))
    stop(error_msg)
  }

  if (!file.exists(output_file)) {
    cat("    Error: Python script failed to create output file:", output_file, "\n")
    # Clean up temp files
    unlink(c(gene_names_file, output_file))
    return(NULL)
  }

  # Read CSV
  tryCatch({
    expr_matrix <- read.csv(output_file, row.names = 1, check.names = FALSE)

    # Check if matrix is empty - this likely means genes weren't found
    if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
      requested_genes <- if (length(gene_names) == 1) gene_names else paste(length(gene_names), "genes")
      error_msg <- paste0(
        "ERROR: Expression matrix is empty - gene(s) not found in RNA data.\n",
        "    Requested: ", requested_genes, "\n",
        "    Period: ", period, "\n",
        "    The pipeline will stop to prevent hanging."
      )
      cat("    ", error_msg, "\n")
      unlink(c(gene_names_file, output_file))
      stop(error_msg)
    }

    # Clean up temp files
    unlink(c(gene_names_file, output_file))

    return(as.matrix(expr_matrix))
  }, error = function(e) {
    cat("    Error reading expression matrix:", e$message, "\n")
    unlink(c(gene_names_file, output_file))
    return(NULL)
  })
}

#' Load metacell assignments for a batch
#' @param batch_name Batch/sample name
#' @return List mapping metacell IDs to cell IDs
load_metacell_assignments <- function(batch_name) {
  mapping_files <- list.files(METACELL_DIR, pattern = paste0("mapping_p[0-9]+_", batch_name, "\\.csv$"),
                              full.names = TRUE)

  if (length(mapping_files) == 0) {
    return(NULL)
  }

  mapping_df <- read.csv(mapping_files[1], stringsAsFactors = FALSE, check.names = FALSE)
  if (!("SEACell" %in% colnames(mapping_df))) {
    return(NULL)
  }

  if ("archr_cell_id_std" %in% colnames(mapping_df)) {
    mapping_df$cell_id_for_match <- as.character(mapping_df$archr_cell_id_std)
  } else if (all(c("library", "cell_id_std") %in% colnames(mapping_df))) {
    mapping_df$cell_id_for_match <- paste0(mapping_df$library, "#", mapping_df$cell_id_std)
  } else if ("index" %in% colnames(mapping_df)) {
    mapping_df$cell_id_for_match <- as.character(mapping_df$index)
  } else {
    return(NULL)
  }

  assignments <- list()
  for (mc_id in unique(mapping_df$SEACell)) {
    cell_ids <- mapping_df$cell_id_for_match[mapping_df$SEACell == mc_id]
    assignments[[mc_id]] <- cell_ids
  }

  return(assignments)
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

load_metacell_assignments_filtered <- function(batch_name, cell_type = NULL) {
  mapping <- load_metacell_assignments(batch_name)
  if (is.null(mapping) || is.null(cell_type) || toupper(cell_type) == "ALL") {
    return(mapping)
  }

  mapping_files <- list.files(METACELL_DIR, pattern = paste0("mapping_p[0-9]+_", batch_name, "\\.csv$"),
                              full.names = TRUE)
  if (length(mapping_files) == 0) {
    return(NULL)
  }
  mapping_df <- read.csv(mapping_files[1], stringsAsFactors = FALSE, check.names = FALSE)

  if ("broad_cell_type" %in% colnames(mapping_df)) {
    mapping_df <- mapping_df[mapping_df$broad_cell_type == cell_type, , drop = FALSE]
  } else if ("majorclass_v1" %in% colnames(mapping_df)) {
    broad <- map_majorclass_to_broad(mapping_df$majorclass_v1)
    mapping_df <- mapping_df[broad == cell_type, , drop = FALSE]
  } else {
    if (cell_type != "ExN") {
      return(NULL)
    }
    # Backward compatibility with older excitatory-only mappings
  }

  if (nrow(mapping_df) == 0) {
    return(NULL)
  }

  if ("archr_cell_id_std" %in% colnames(mapping_df)) {
    mapping_df$cell_id_for_match <- as.character(mapping_df$archr_cell_id_std)
  } else if (all(c("library", "cell_id_std") %in% colnames(mapping_df))) {
    mapping_df$cell_id_for_match <- paste0(mapping_df$library, "#", mapping_df$cell_id_std)
  } else if ("index" %in% colnames(mapping_df)) {
    mapping_df$cell_id_for_match <- as.character(mapping_df$index)
  } else {
    return(NULL)
  }

  assignments <- list()
  for (mc_id in unique(mapping_df$SEACell)) {
    cell_ids <- mapping_df$cell_id_for_match[mapping_df$SEACell == mc_id]
    assignments[[mc_id]] <- cell_ids
  }
  assignments
}

#' Aggregate RNA expression to metacell level (same logic as batch version)
aggregate_rna_to_metacells <- function(expr_matrix, metacell_assignments) {
  expr_metacells <- rownames(expr_matrix)

  n_metacells <- length(metacell_assignments)
  n_genes <- ncol(expr_matrix)
  metacell_expr <- matrix(0, nrow = n_metacells, ncol = n_genes)
  rownames(metacell_expr) <- names(metacell_assignments)
  colnames(metacell_expr) <- colnames(expr_matrix)

  pb <- txtProgressBar(min = 0, max = n_metacells, style = 3, char = "=")
  matched_count <- 0

  for (i in seq_along(metacell_assignments)) {
    mc_id <- names(metacell_assignments)[i]

    mc_id_simple <- mc_id
    if (grepl("_seacell_", mc_id)) {
      mc_id_simple <- strsplit(mc_id, "_seacell_")[[1]][2]
    } else if (grepl("_", mc_id)) {
      parts <- strsplit(mc_id, "_")[[1]]
      if (length(parts) > 0 && grepl("SEACell", parts[length(parts)])) {
        mc_id_simple <- parts[length(parts)]
      }
    }

    matched_idx <- NULL
    if (mc_id %in% expr_metacells) {
      matched_idx <- which(expr_metacells == mc_id)[1]
    } else if (mc_id_simple %in% expr_metacells) {
      matched_idx <- which(expr_metacells == mc_id_simple)[1]
    } else {
      pattern <- paste0("^", gsub("-", "\\-", mc_id_simple, fixed = TRUE), "-0$")
      partial_match <- grep(pattern, expr_metacells, value = FALSE)
      if (length(partial_match) > 0) {
        matched_idx <- partial_match[1]
      } else {
        prefix_match <- which(startsWith(expr_metacells, paste0(mc_id_simple, "-")))
        if (length(prefix_match) > 0) {
          matched_idx <- prefix_match[1]
        } else {
          expr_ids_trimmed <- gsub("-0$", "", expr_metacells)
          trimmed_match <- which(expr_ids_trimmed == mc_id_simple)
          if (length(trimmed_match) > 0) {
            matched_idx <- trimmed_match[1]
          }
        }
      }
    }

    if (!is.null(matched_idx)) {
      matched_count <- matched_count + 1
      expr_row <- expr_matrix[matched_idx, , drop = FALSE]
      if (inherits(expr_row, "sparseMatrix")) {
        metacell_expr[i, ] <- as.numeric(as.matrix(expr_row))
      } else {
        metacell_expr[i, ] <- as.numeric(expr_row)
      }
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(metacell_expr)
}

#' Load ArchRProject for a specific batch from Arrow file
load_batch_archr_from_arrow <- function(batch_name) {
  arrow_file_pattern <- paste0(".*ArrowFiles-", batch_name, "\\.arrow$")
  arrow_files <- list.files(ATAC_DIR, pattern = arrow_file_pattern, full.names = TRUE)

  if (length(arrow_files) == 0) {
    stop("ERROR: Arrow file not found for batch ", batch_name, " in ", ATAC_DIR)
  }

  arrow_file <- arrow_files[1]
  temp_proj_dir <- file.path(BASE_DIR, "temp_archr_projects", batch_name)
  if (dir.exists(temp_proj_dir)) {
    unlink(temp_proj_dir, recursive = TRUE)
  }
  dir.create(temp_proj_dir, recursive = TRUE)

  proj_batch <- ArchRProject(
    ArrowFiles = arrow_file,
    outputDirectory = temp_proj_dir,
    copyArrows = FALSE
  )

  peak_set <- get_main_peak_set_cached()
  if (!is.null(peak_set) && length(peak_set) > 0) {
      proj_batch <- addPeakSet(ArchRProj = proj_batch, peakSet = peak_set, force = TRUE)

      need_peak_matrix <- TRUE
      available_matrices <- getAvailableMatrices(proj_batch)
      if ("PeakMatrix" %in% available_matrices) {
        peak_se <- tryCatch(
          getMatrixFromProject(ArchRProj = proj_batch, useMatrix = "PeakMatrix"),
          error = function(e) NULL
        )
        if (!is.null(peak_se)) {
          pm <- assay(peak_se)
          need_peak_matrix <- (nrow(pm) == 0 || ncol(pm) == 0)
        }
      }

      if (need_peak_matrix) {
        proj_batch <- addPeakMatrix(proj_batch)
      }
  }

  return(proj_batch)
}

#' Aggregate ATAC data to metacell level (same logic as batch version)
aggregate_atac_to_metacells <- function(proj, batch_name, peak_coords, metacell_assignments) {
  peak_set <- getPeakSet(proj)
  if (is.null(peak_set) || length(peak_set) == 0) {
    if (file.exists(MAIN_ARCHR_PROJ_FILE)) {
      main_proj <- readRDS(MAIN_ARCHR_PROJ_FILE)
      peak_set <- getPeakSet(main_proj)
      if (!is.null(peak_set) && length(peak_set) > 0) {
        proj <- addPeakSet(ArchRProj = proj, peakSet = peak_set, force = TRUE)
      }
    }
    if (is.null(peak_set) || length(peak_set) == 0) {
      stop("ERROR: Peak set not found in ArchRProject")
    }
  }

  peak_granges <- do.call(c, lapply(peak_coords, function(coord) {
    parts <- strsplit(coord, "-")[[1]]
    GRanges(seqnames = parts[1],
            ranges = IRanges(start = as.numeric(parts[2]), end = as.numeric(parts[3])))
  }))

  overlaps <- findOverlaps(peak_set, peak_granges)
  matched_peaks <- peak_set[queryHits(overlaps)]
  matched_coords <- paste(seqnames(matched_peaks), start(matched_peaks), end(matched_peaks), sep = "-")

  if (length(matched_coords) == 0) {
    stop("ERROR: No peaks matched between motif results and ArchRProject")
  }

  batch_cells <- proj$cellNames
  if (length(batch_cells) == 0) {
    stop("ERROR: No cells found in ArchRProject")
  }

  peak_matrix_se <- getMatrixFromProject(ArchRProj = proj, useMatrix = "PeakMatrix")
  peak_matrix <- assay(peak_matrix_se)

  # Get peak coordinates from peak matrix rownames
  # ArchR peak matrix rownames might use dots (chr12.24469992.24470492) instead of dashes
  peak_matrix_rownames <- rownames(peak_matrix)

  # Normalize peak matrix rownames to use dashes (to match motif matrix format)
  # Convert dots to dashes: chr12.24469992.24470492 -> chr12-24469992-24470492
  peak_matrix_coords_normalized <- gsub("\\.", "-", peak_matrix_rownames)

  # Create peak_set_coords with dashes (for matching)
  peak_set_coords <- paste(seqnames(peak_set), start(peak_set), end(peak_set), sep = "-")

  # Try to match using normalized coordinates
  # First try matching against normalized peak matrix rownames
  matched_peak_indices <- match(matched_coords, peak_matrix_coords_normalized)
  matched_peak_indices <- matched_peak_indices[!is.na(matched_peak_indices)]

  if (length(matched_peak_indices) == 0) {
    # Try matching against peak_set_coords (which already uses dashes)
    matched_peak_indices <- match(matched_coords, peak_set_coords)
    matched_peak_indices <- matched_peak_indices[!is.na(matched_peak_indices)]
  }

  if (length(matched_peak_indices) == 0) {
    # Debug: show sample names for troubleshooting
    cat("  Debug: Sample matched_coords:", paste(head(matched_coords, 3), collapse = ", "), "\n")
    cat("  Debug: Sample peak_matrix_coords_normalized:", paste(head(peak_matrix_coords_normalized, 3), collapse = ", "), "\n")
    cat("  Debug: Sample peak_set_coords:", paste(head(peak_set_coords, 3), collapse = ", "), "\n")
    stop("ERROR: Could not match peak coordinates to peak matrix rows")
  }

  peak_subset <- peak_matrix[matched_peak_indices, , drop = FALSE]
  # Ensure rownames use dashes (consistent with motif matrix format: chr-start-end)
  # Use matched_coords which already has the correct format with dashes
  rownames(peak_subset) <- matched_coords

  n_metacells <- length(metacell_assignments)
  n_peaks <- length(matched_coords)
  atac_metacell_matrix <- matrix(0, nrow = n_metacells, ncol = n_peaks)
  rownames(atac_metacell_matrix) <- names(metacell_assignments)
  colnames(atac_metacell_matrix) <- matched_coords

  pb <- txtProgressBar(min = 0, max = length(metacell_assignments), style = 3, char = "=")

  for (i in seq_along(metacell_assignments)) {
    mc_id <- names(metacell_assignments)[i]
    cell_ids <- metacell_assignments[[mc_id]]

    peak_cols <- colnames(peak_subset)
    map_ids <- unique(canonical_cell_id(cell_ids, sample_hint = batch_name))
    peak_ids <- canonical_cell_id(peak_cols, sample_hint = batch_name)
    matched_cells <- peak_cols[peak_ids %in% map_ids]

    if (length(matched_cells) > 0) {
      if (length(matched_cells) == 1) {
        atac_metacell_matrix[i, ] <- as.numeric(peak_subset[, matched_cells])
      } else {
        atac_metacell_matrix[i, ] <- Matrix::rowSums(peak_subset[, matched_cells, drop = FALSE])
      }
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  return(atac_metacell_matrix)
}

# =============================================================================
# BATCH CORRECTION FUNCTIONS
# =============================================================================

#' Apply Seurat RPCA batch correction for RNA data
#' @param rna_data RNA expression matrix (samples × genes)
#' @param batch_info Batch labels for each sample
#' @return Batch-corrected RNA matrix
apply_rna_batch_correction <- function(rna_data, batch_info) {
  cat("  Using Seurat RPCA for RNA batch correction\n")

  # Ensure data is a numeric matrix (not data frame)
  if (inherits(rna_data, "data.frame")) {
    rna_data <- as.matrix(rna_data)
  }
  if (is.null(rownames(rna_data))) {
    rownames(rna_data) <- paste0("Cell_", seq_len(nrow(rna_data)))
  }

  # Convert sparse matrix to regular matrix if needed
  if (inherits(rna_data, "sparseMatrix")) {
    rna_data <- as.matrix(rna_data)
  }

  # Ensure it's numeric and finite
  storage.mode(rna_data) <- "numeric"
  rna_data[is.na(rna_data)] <- 0
  rna_data[is.infinite(rna_data)] <- 0

  # Transpose: Seurat expects genes × cells format
  rna_data_t <- t(rna_data)

  # Ensure row and column names are set (required by Seurat)
  if (is.null(rownames(rna_data_t))) {
    rownames(rna_data_t) <- paste0("Gene_", seq_len(nrow(rna_data_t)))
  }
  if (is.null(colnames(rna_data_t))) {
    colnames(rna_data_t) <- paste0("Cell_", seq_len(ncol(rna_data_t)))
  }

  # Convert to sparse matrix explicitly for Seurat v5 compatibility
  # Seurat v5 prefers sparse matrices and has issues with dense matrices
  if (!inherits(rna_data_t, "sparseMatrix")) {
    # Convert to sparse matrix using Matrix package
    rna_data_t <- Matrix::Matrix(rna_data_t, sparse = TRUE)
  }

  # Ensure the matrix is properly formatted (2D)
  if (length(dim(rna_data_t)) != 2) {
    stop("ERROR: Data must be a 2-dimensional matrix")
  }

  # Ensure it's a dgCMatrix (the format Seurat expects)
  if (!inherits(rna_data_t, "dgCMatrix")) {
    rna_data_t <- as(rna_data_t, "CsparseMatrix")
  }

  cell_names <- colnames(rna_data_t)
  batch_info <- as.character(batch_info)
  names(batch_info) <- cell_names

  # Split by batch and run RPCA integration
  cells_by_batch <- split(cell_names, batch_info[cell_names])
  obj_list <- lapply(cells_by_batch, function(cells) {
    obj <- suppressWarnings(Seurat::CreateSeuratObject(counts = rna_data_t[, cells, drop = FALSE]))
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
    return(obj)
  })

  n_batches <- length(obj_list)
  if (n_batches <= 1) {
    corrected_data <- t(Seurat::GetAssayData(obj_list[[1]], layer = "data"))
  } else {
    min_batch_cells <- min(sapply(obj_list, ncol))
    n_features <- nrow(rna_data_t)
    if (min(n_features, min_batch_cells) <= 2) {
      stop("ERROR: Insufficient dimensions for RPCA (features=", n_features,
           ", min batch cells=", min_batch_cells, ")")
    }

    features <- Seurat::SelectIntegrationFeatures(object.list = obj_list)
    all_input_features <- rownames(rna_data_t)
    common_features <- Reduce(intersect, lapply(obj_list, rownames))
    features_to_integrate <- intersect(all_input_features, common_features)
    if (length(features_to_integrate) == 0) {
      stop("ERROR: No common features available for RPCA integration")
    }
    cat("  RPCA features: anchors =", length(features),
        ", integrate =", length(features_to_integrate), "\n")
    obj_list <- lapply(obj_list, function(obj) {
      obj <- Seurat::ScaleData(obj, features = features, verbose = FALSE)
      obj <- Seurat::RunPCA(obj, features = features, verbose = FALSE)
      obj
    })

    max_dims <- min(30, min(sapply(obj_list, function(x) ncol(x))) - 1)
    if (!is.finite(max_dims) || max_dims < 2) {
      max_dims <- 2
    }
    anchors <- Seurat::FindIntegrationAnchors(
      object.list = obj_list,
      anchor.features = features,
      reduction = "rpca",
      dims = 1:max_dims
    )
    k_weight <- max(10, min(50, floor(min_batch_cells / 2)))
    cat("  RPCA integration params: dims =", max_dims, ", k.weight =", k_weight, "\n")
    integrated <- Seurat::IntegrateData(
      anchorset = anchors,
      dims = 1:max_dims,
      k.weight = k_weight,
      features.to.integrate = features_to_integrate
    )
    corrected_data <- t(Seurat::GetAssayData(integrated, assay = "integrated", layer = "data"))
  }

  # Convert back to regular matrix if needed
  if (inherits(corrected_data, "sparseMatrix")) {
    corrected_data <- as.matrix(corrected_data)
  }

  # Ensure original sample order
  if (!is.null(rownames(corrected_data)) && !is.null(rownames(rna_data))) {
    corrected_data <- corrected_data[rownames(rna_data), , drop = FALSE]
  }

  return(corrected_data)
}

#' Apply ATAC batch correction using TF-IDF + LSI embeddings
#' @param atac_data ATAC accessibility matrix (samples × peaks)
#' @param batch_info Batch labels for each sample
#' @return Batch-corrected ATAC matrix
apply_atac_batch_correction <- function(atac_data, batch_info) {
  cat("  Using LSI embedding correction for ATAC batch correction\n")

  atac_data <- as.matrix(atac_data)
  atac_data[is.na(atac_data)] <- 0
  atac_data[is.infinite(atac_data)] <- 0
  atac_data[atac_data < 0] <- 0

  # TF-IDF transform (cell x peak)
  rs <- rowSums(atac_data)
  rs[rs == 0] <- 1
  tf <- atac_data / rs

  df <- colSums(atac_data > 0)
  n_cells <- nrow(atac_data)
  idf <- log(1 + n_cells / (1 + df))
  tfidf <- sweep(tf, 2, idf, FUN = "*")

  k <- min(30, nrow(tfidf) - 1, ncol(tfidf) - 1)
  if (k < 2) {
    return(atac_data)
  }

  sv <- irlba::irlba(tfidf, nv = k, nu = k)
  emb <- sv$u %*% diag(sv$d)  # cell x k

  # Batch centering in LSI space
  batch_info <- as.character(batch_info)
  global_mean <- colMeans(emb)
  emb_corr <- emb
  for (b in unique(batch_info)) {
    idx <- which(batch_info == b)
    if (length(idx) == 0) next
    bmean <- colMeans(emb[idx, , drop = FALSE])
    emb_corr[idx, ] <- sweep(emb[idx, , drop = FALSE], 2, bmean, "-")
    emb_corr[idx, ] <- sweep(emb_corr[idx, , drop = FALSE], 2, global_mean, "+")
  }

  # Reconstruct corrected accessibility-like matrix for downstream modeling
  reconstructed <- emb_corr %*% t(sv$v)
  reconstructed[!is.finite(reconstructed)] <- 0
  reconstructed[reconstructed < 0] <- 0
  rownames(reconstructed) <- rownames(atac_data)
  colnames(reconstructed) <- colnames(atac_data)

  return(reconstructed)
}

# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

#' Run period-specific data retrieval pipeline
#' @param target_gene Target gene symbol
#' @param period Period number
#' @param output_dir Output directory
#' @param apply_batch_correction Whether to apply batch correction (RNA: Seurat RPCA, ATAC: LSI) (default: TRUE)
#' @return List with corrected X, Y, Z matrices
run_period_data_retrieval <- function(target_gene, period, output_dir = OUTPUT_DIR,
                                     apply_batch_correction = TRUE,
                                     cell_type = NULL,
                                     motif_results_dir = PERIOD_MOTIF_RESULTS_DIR,
                                     n_cores = 10) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PERIOD-SPECIFIC DATA RETRIEVAL PIPELINE\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Target Gene:", target_gene, "\n")
  cat("Period:", period, "\n")
  cat("Cell Type:", ifelse(is.null(cell_type), "ALL", cell_type), "\n")
  cat("Apply Batch Correction:", apply_batch_correction, "\n")
  cat("Requested cores:", n_cores, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  addArchRThreads(threads = min(4, max(1, as.integer(n_cores))))

  # Create output directory
  cell_type_suffix <- if (is.null(cell_type) || toupper(cell_type) == "ALL") "all_cells" else paste0("celltype_", cell_type)
  period_output_dir <- file.path(output_dir, paste0("period_", period), cell_type_suffix, target_gene)
  if (!dir.exists(period_output_dir)) {
    dir.create(period_output_dir, recursive = TRUE)
  }

  # Step 1: Check period motif matching results
  cat("Step 1: Checking period motif matching results...\n")
  motif_check <- check_period_motif_results(target_gene, period, motif_results_dir = motif_results_dir)

  if (!motif_check$exists) {
    error_msg <- paste("ERROR:", motif_check$message)
    cat(error_msg, "\n")
    stop(error_msg)
  }

  cat("✓", motif_check$message, "\n")

  # Step 2: Load motif matrix
  cat("\nStep 2: Loading motif matrix...\n")
  motif_matrix <- read.csv(motif_check$motif_matrix_file, row.names = 1, check.names = FALSE)
  cat("✓ Loaded motif matrix:", dim(motif_matrix), "\n")

  tf_names <- rownames(motif_matrix)
  peak_coords <- colnames(motif_matrix)

  # Step 3: Get batches in period
  cat("\nStep 3: Identifying batches in period...\n")
  batches <- get_batches_for_period(period)

  if (length(batches) == 0) {
    stop("ERROR: No batches found for period ", period)
  }

  cat("✓ Found", length(batches), "batches:", paste(batches, collapse = ", "), "\n")

  # Step 4: Load RNA expression once for the entire period (shared across batches)
  cat("\nStep 4: Loading RNA expression for period (shared across all batches)...\n")
  X_matrix_cells <- load_rna_expression_for_period(period, tf_names)

  if (is.null(X_matrix_cells)) {
    stop("ERROR: Could not load TF expression (X matrix) for period ", period,
         "\nThis affects all batches. Please check:\n",
         "  1. h5ad file exists: ", file.path(METACELL_DIR, paste0("metacells_period_", period, ".h5ad")), "\n",
         "  2. Python script exists: ", PYTHON_EXTRACT_RNA_SCRIPT, "\n",
         "  3. TF gene names are correct")
  }

  # Load target gene expression - stop immediately if not found
  Y_matrix_cells <- load_rna_expression_for_period(period, target_gene)

  if (is.null(Y_matrix_cells)) {
    # Provide helpful error message with common alternative symbols
    alternative_hints <- ""
    if (target_gene == "SUV420H1") {
      alternative_hints <- "\n  Note: SUV420H1 may be stored as 'KMT5B' in the expression data."
    } else if (target_gene == "SUV420H2") {
      alternative_hints <- "\n  Note: SUV420H2 may be stored as 'KMT5C' in the expression data."
    }

    error_msg <- paste0(
      "ERROR: Target gene '", target_gene, "' not found in RNA expression data for period ", period, ".",
      alternative_hints,
      "\n  Please check the gene symbol used in the expression data (h5ad file).",
      "\n  The pipeline will stop to prevent hanging."
    )
    cat(error_msg, "\n")
    stop(error_msg)
  }

  # Remove target gene from X before any downstream aggregation/correction/saving.
  x_cols_upper <- toupper(colnames(X_matrix_cells))
  target_upper <- toupper(target_gene)
  target_in_x <- which(x_cols_upper == target_upper)
  if (length(target_in_x) > 0) {
    removed_cols <- colnames(X_matrix_cells)[target_in_x]
    X_matrix_cells <- X_matrix_cells[, -target_in_x, drop = FALSE]
    cat("✓ Removed target gene from X prior to data retrieval:", paste(removed_cols, collapse = ", "), "\n")
  } else {
    cat("✓ Target gene not present in X prior to data retrieval (no removal needed)\n")
  }
  if (ncol(X_matrix_cells) == 0) {
    stop("ERROR: X matrix has no TF predictors after removing target gene")
  }

  cat("✓ Loaded RNA expression:", nrow(X_matrix_cells), "metacells ×", ncol(X_matrix_cells), "TFs\n")
  cat("✓ Loaded target gene expression:", nrow(Y_matrix_cells), "metacells\n")

  # Step 5: Retrieve data for each batch
  cat("\nStep 5: Retrieving data for each batch...\n")

  batch_cache_dir <- file.path(period_output_dir, "batch_cache")
  if (!dir.exists(batch_cache_dir)) {
    dir.create(batch_cache_dir, recursive = TRUE)
  }
  tf_sig <- compact_signature(tf_names)
  peak_sig <- compact_signature(peak_coords)

  X_batch_list <- list()
  Y_batch_list <- list()
  Z_batch_list <- list()
  batch_labels <- character()

  for (batch_idx in seq_along(batches)) {
    batch_name <- batches[batch_idx]
    cat("\n  Processing batch", batch_idx, "of", length(batches), ":", batch_name, "\n")

    cache_prefix <- file.path(batch_cache_dir, batch_name)
    cache_x <- paste0(cache_prefix, "_X.rds")
    cache_y <- paste0(cache_prefix, "_Y.rds")
    cache_z <- paste0(cache_prefix, "_Z.rds")
    cache_meta <- paste0(cache_prefix, "_meta.rds")
    can_use_cache <- FALSE
    if (file.exists(cache_x) && file.exists(cache_y) && file.exists(cache_z) && file.exists(cache_meta)) {
      meta <- tryCatch(readRDS(cache_meta), error = function(e) NULL)
      can_use_cache <- !is.null(meta) &&
        identical(meta$target_gene, target_gene) &&
        identical(meta$period, as.integer(period)) &&
        identical(meta$cell_type, ifelse(is.null(cell_type), "ALL", cell_type)) &&
        identical(meta$tf_sig, tf_sig) &&
        identical(meta$peak_sig, peak_sig) &&
        isTRUE(meta$target_removed_from_x)
    }
    if (can_use_cache) {
      X_cached <- tryCatch(readRDS(cache_x), error = function(e) NULL)
      Y_cached <- tryCatch(readRDS(cache_y), error = function(e) NULL)
      Z_cached <- tryCatch(readRDS(cache_z), error = function(e) NULL)
      if (!is.null(X_cached) && !is.null(Y_cached) && !is.null(Z_cached) &&
          nrow(X_cached) > 0 && nrow(Y_cached) > 0 && nrow(Z_cached) > 0) {
        X_batch_list[[batch_name]] <- X_cached
        Y_batch_list[[batch_name]] <- Y_cached
        Z_batch_list[[batch_name]] <- Z_cached
        batch_labels <- c(batch_labels, rep(batch_name, nrow(X_cached)))
        cat("    ✓ Loaded cached batch matrices\n")
        next
      }
    }

    # Load metacell assignments
    metacell_assignments <- load_metacell_assignments_filtered(batch_name, cell_type = cell_type)
    if (is.null(metacell_assignments) || length(metacell_assignments) == 0) {
      cat("  Warning: Skipping batch", batch_name, "- no metacell assignments for selected cell type\n")
      next
    }

    # Aggregate to metacell level
    X_batch <- aggregate_rna_to_metacells(X_matrix_cells, metacell_assignments)
    Y_batch <- aggregate_rna_to_metacells(Y_matrix_cells, metacell_assignments)

    # Add batch prefix to metacell IDs to avoid conflicts
    rownames(X_batch) <- paste0(batch_name, "_", rownames(X_batch))
    rownames(Y_batch) <- paste0(batch_name, "_", rownames(Y_batch))

    # Load ATAC data
    cat("    Loading ATAC data...\n")
    tryCatch({
      proj_batch <- load_batch_archr_from_arrow(batch_name)
      Z_batch <- aggregate_atac_to_metacells(proj_batch, batch_name, peak_coords, metacell_assignments)
      rownames(Z_batch) <- paste0(batch_name, "_", rownames(Z_batch))
      # Normalize column names to use dashes (consistent with motif matrix format)
      colnames(Z_batch) <- gsub("\\.", "-", colnames(Z_batch))
    }, error = function(e) {
      cat("  Warning: Skipping batch", batch_name, "- ATAC data error:", e$message, "\n")
      Z_batch <<- NULL
    })

    if (is.null(Z_batch)) {
      next
    }

    # Align matrices by common metacells
    common_metacells <- intersect(intersect(rownames(X_batch), rownames(Y_batch)), rownames(Z_batch))

    if (length(common_metacells) == 0) {
      cat("  Warning: Skipping batch", batch_name, "- no common metacells\n")
      next
    }

    X_batch <- X_batch[common_metacells, , drop = FALSE]
    Y_batch <- Y_batch[common_metacells, , drop = FALSE]
    Z_batch <- Z_batch[common_metacells, , drop = FALSE]

    # Store matrices and batch labels
    X_batch_list[[batch_name]] <- X_batch
    Y_batch_list[[batch_name]] <- Y_batch
    Z_batch_list[[batch_name]] <- Z_batch
    batch_labels <- c(batch_labels, rep(batch_name, nrow(X_batch)))

    saveRDS(X_batch, cache_x)
    saveRDS(Y_batch, cache_y)
    saveRDS(Z_batch, cache_z)
    saveRDS(list(
      target_gene = target_gene,
      period = as.integer(period),
      cell_type = ifelse(is.null(cell_type), "ALL", cell_type),
      tf_sig = tf_sig,
      peak_sig = peak_sig,
      n_metacells = nrow(X_batch),
      target_removed_from_x = TRUE
    ), cache_meta)

    cat("    ✓ Retrieved", nrow(X_batch), "metacells from batch", batch_name, "\n")
  }

  if (length(X_batch_list) == 0) {
    stop("ERROR: No batches successfully processed")
  }

  cat("\n✓ Retrieved data from", length(X_batch_list), "batches\n")

  # Step 5: Combine matrices from all batches
  cat("\nStep 5: Combining matrices from all batches...\n")

  X_combined <- do.call(rbind, X_batch_list)
  Y_combined <- do.call(rbind, Y_batch_list)
  Z_combined <- do.call(rbind, Z_batch_list)

  # Normalize Z matrix column names to use dashes (consistent with motif matrix format)
  # ArchR might use dots (chr12.24469992.24470492), convert to dashes (chr12-24469992-24470492)
  colnames(Z_combined) <- gsub("\\.", "-", colnames(Z_combined))

  cat("✓ Combined matrices:\n")
  cat("  X:", dim(X_combined), "\n")
  cat("  Y:", dim(Y_combined), "\n")
  cat("  Z:", dim(Z_combined), "\n")
  cat("  Z column names normalized to use dashes (matching motif matrix format)\n")

  # Step 6: Batch correction
  cat("\nStep 6: Batch correction...\n")

  # Check if batch correction is needed (requires multiple batches)
  n_unique_batches <- length(unique(batch_labels))
  cat("  Number of unique batches:", n_unique_batches, "\n")

  if (!apply_batch_correction) {
    cat("  Batch correction disabled - using uncorrected data\n")
    X_corrected <- X_combined
    Y_corrected <- Y_combined
    Z_corrected <- Z_combined
    batch_correction_applied <- FALSE
  } else if (n_unique_batches <= 1) {
    cat("  Only one batch found - skipping batch correction\n")
    X_corrected <- X_combined
    Y_corrected <- Y_combined
    Z_corrected <- Z_combined
    batch_correction_applied <- FALSE
  } else {
    cat("  Multiple batches found - applying batch correction (RNA: RPCA, ATAC: LSI)...\n")

    # Combine X and Y for batch correction (Y is one-dimensional, so combine them)
    cat("  Combining X and Y matrices for joint batch correction...\n")
    # Use a temporary unique target column name to avoid feature-name collisions
    # when target gene is also present in X (e.g., FEZF2 in TF list).
    # Use alphanumeric-only feature name to avoid Seurat feature-name sanitization.
    temp_target_col <- "TARGETGENEFORRPCA"
    while (temp_target_col %in% colnames(X_combined)) {
      temp_target_col <- paste0(temp_target_col, "_X")
    }
    Y_for_correction <- Y_combined
    colnames(Y_for_correction) <- temp_target_col
    XY_combined <- cbind(X_combined, Y_for_correction)

    # Seurat requires unique feature names. Keep a lookup to restore original names later.
    x_lookup <- colnames(X_combined)
    target_lookup <- temp_target_col
    if (anyDuplicated(colnames(XY_combined)) > 0) {
      dup_names <- unique(colnames(XY_combined)[duplicated(colnames(XY_combined))])
      cat("  Warning: Duplicate feature names detected before RPCA; applying make.unique().\n")
      cat("    Example duplicates:", paste(head(dup_names, 10), collapse = ", "), "\n")
      colnames(XY_combined) <- make.unique(colnames(XY_combined), sep = "__dup")
      x_lookup <- colnames(XY_combined)[seq_len(ncol(X_combined))]
      target_lookup <- colnames(XY_combined)[ncol(XY_combined) + 1]
    }
    cat("  Combined XY matrix:", dim(XY_combined), "\n")

    # Apply batch correction to combined XY matrix
    # Fallback: for very small batches/samples, RPCA can fail due to irlba dimension limits.
    XY_corrected <- tryCatch(
      apply_rna_batch_correction(XY_combined, batch_labels),
      error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("max\\(nu, nv\\) must be strictly less than min\\(nrow\\(A\\), ncol\\(A\\)\\)", msg) ||
            grepl("Insufficient dimensions for RPCA", msg)) {
          cat("  Warning: RPCA skipped due to small sample size; using uncorrected RNA for this job.\n")
          cat("    RPCA error:", msg, "\n")
          return(XY_combined)
        }
        stop(e)
      }
    )

    # Split back into X and Y by feature name (RPCA can reorder columns)
    xy_cols <- colnames(XY_corrected)
    if (is.null(xy_cols)) {
      stop("ERROR: Batch-corrected XY matrix has no column names; cannot map X/Y by feature name")
    }

    missing_x <- setdiff(x_lookup, xy_cols)
    if (length(missing_x) > 0) {
      stop("ERROR: ", length(missing_x), " TF columns missing after RNA batch correction (e.g., ",
           paste(head(missing_x, 5), collapse = ", "), ")")
    }

    if (!(target_lookup %in% xy_cols)) {
      # Seurat may sanitize feature names (e.g., underscores -> dashes).
      alt_target_lookup <- gsub("_", "-", target_lookup, fixed = TRUE)
      if (alt_target_lookup %in% xy_cols) {
        cat("  Warning: target feature name was sanitized by Seurat; using:", alt_target_lookup, "\n")
        target_lookup <- alt_target_lookup
      } else {
        stop("ERROR: Target feature column '", target_lookup, "' not found in batch-corrected XY matrix")
      }
    }

    X_corrected <- XY_corrected[, x_lookup, drop = FALSE]
    colnames(X_corrected) <- colnames(X_combined)
    Y_corrected <- XY_corrected[, target_lookup, drop = FALSE]
    colnames(Y_corrected) <- target_gene

    cat("  ✓ Split corrected XY back into X and Y\n")
    cat("    X corrected:", dim(X_corrected), "\n")
    cat("    Y corrected:", dim(Y_corrected), "\n")

    # Apply batch correction to Z (ATAC)
    Z_corrected <- apply_atac_batch_correction(Z_combined, batch_labels)

    cat("✓ Applied batch correction to all matrices\n")
    batch_correction_applied <- TRUE
  }

  # Step 7: Save results
  cat("\nStep 7: Saving results...\n")

  write.csv(X_corrected, file.path(period_output_dir, paste0(target_gene, "_X_matrix.csv")))
  write.csv(Y_corrected, file.path(period_output_dir, paste0(target_gene, "_Y_matrix.csv")))
  write.csv(Z_corrected, file.path(period_output_dir, paste0(target_gene, "_Z_matrix.csv")))
  write.csv(motif_matrix, file.path(period_output_dir, paste0(target_gene, "_motif_matrix.csv")))

  # Save batch labels (always save, regardless of batch correction status)
  batch_info_df <- data.frame(
    metacell_id = rownames(X_corrected),
    batch = batch_labels,
    stringsAsFactors = FALSE
  )
  write.csv(batch_info_df, file.path(period_output_dir, paste0(target_gene, "_batch_info.csv")), row.names = FALSE)
  cat("✓ Saved batch information to:", file.path(period_output_dir, paste0(target_gene, "_batch_info.csv")), "\n")

  # Save summary
  summary_file <- file.path(period_output_dir, paste0(target_gene, "_data_retrieval_summary.txt"))

  # Batch correction status is already determined above

  summary_lines <- c(
    paste("Period-Specific Data Retrieval Pipeline Summary"),
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Period:", period),
    paste("Cell Type:", ifelse(is.null(cell_type), "ALL", cell_type)),
    paste("Batches Processed:", length(X_batch_list)),
    paste("Batches:", paste(names(X_batch_list), collapse = ", ")),
    paste("Unique Batches:", n_unique_batches),
    paste(""),
    paste("Matrix Dimensions:"),
    paste("  X (TFs):", nrow(X_corrected), "metacells ×", ncol(X_corrected), "TFs"),
    paste("  Y (Target):", nrow(Y_corrected), "metacells ×", ncol(Y_corrected), "target gene"),
    paste("  Z (ATAC):", nrow(Z_corrected), "metacells ×", ncol(Z_corrected), "peaks"),
    paste(""),
    paste("Batch Correction:"),
    if (batch_correction_applied) {
      c(paste("  RNA (X, Y): Seurat RPCA (combined correction)"),
        paste("  ATAC (Z): LSI embedding correction"))
    } else if (!apply_batch_correction) {
      c(paste("  RNA (X, Y): None (batch correction disabled)"),
        paste("  ATAC (Z): None (batch correction disabled)"))
    } else {
      c(paste("  RNA (X, Y): None (only one batch)"),
        paste("  ATAC (Z): None (only one batch)"))
    },
    paste(""),
    paste("Output Directory:", period_output_dir)
  )
  writeLines(summary_lines, summary_file)
  cat("✓ Saved summary to:", summary_file, "\n")

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PERIOD-SPECIFIC DATA RETRIEVAL COMPLETED SUCCESSFULLY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  return(list(
    X_matrix = X_corrected,
    Y_matrix = Y_corrected,
    Z_matrix = Z_corrected,
    motif_matrix = motif_matrix,
    batch_info = batch_labels
  ))
}
