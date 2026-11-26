#!/usr/bin/env Rscript

# =============================================================================
# MOTIF MATCHING HELPERS
# =============================================================================
# This script performs motif matching for a target gene at the period level.
# Since motif matching results are identical across batches (same peak set,
# same sequences), this runs motif matching once per period using the main
# ArchRProject and identifies all batches belonging to the period.
#
# Internal module shared by the phase and locus entry points.

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(ArchR)
  library(biomaRt)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2024)
  library(TFBSTools)
  library(Biostrings)
})

# Initialize ArchR only if not being sourced
if (!exists(".SOURCING_PIPELINES")) {
  addArchRThreads(threads = 4)
  addArchRGenome("hg38")
}

# =============================================================================
# MOTIF PIPELINE FUNCTIONS (same as batch version)
# =============================================================================

#' Write MEME motif format from PFM
#' @param pfm Position Frequency Matrix
#' @param file Output file path
write_meme_motif <- function(pfm, file) {
  con <- file(file, "wt")
  on.exit(close(con))

  writeLines("MEME version 4", con)
  writeLines("", con)
  writeLines("ALPHABET= A C G T", con)
  writeLines("", con)
  writeLines("strands: + -", con)
  writeLines("", con)
  writeLines("Background letter frequencies:", con)
  writeLines("A 0.25 C 0.25 G 0.25 T 0.25", con)
  writeLines("", con)

  name <- ID(pfm)
  matrix <- Matrix(pfm)
  matrix <- matrix / colSums(matrix)

  writeLines(sprintf("MOTIF %s", name), con)
  writeLines(sprintf("letter-probability matrix: alength= 4 w= %d nsites= 20 E= 0", ncol(matrix)), con)

  for (i in 1:ncol(matrix)) {
    probs <- sprintf("%.6f", matrix[, i])
    writeLines(paste(probs, collapse = " "), con)
  }
}

get_fasta_get_markov_bin <- function() {
  candidates <- c(
    Sys.getenv("FASTA_GET_MARKOV_BIN", unset = ""),
    Sys.which("fasta-get-markov")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (p in candidates) {
    if (file.exists(p)) return(p)
  }
  ""
}

get_fimo_bin <- function() {
  candidates <- c(
    Sys.getenv("FIMO_BIN", unset = ""),
    Sys.which("fimo")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  for (path in candidates) {
    if (file.exists(path)) return(path)
  }
  stop("FIMO was not found. Add it to PATH or set FIMO_BIN to the executable path.")
}

.TIMEOUT_CMD_CACHE <- NULL

get_timeout_command <- function() {
  if (!is.null(.TIMEOUT_CMD_CACHE)) {
    return(.TIMEOUT_CMD_CACHE)
  }
  timeout_check <- tryCatch(system2("which", "timeout", stdout = TRUE, stderr = TRUE), error = function(e) "")
  if (length(timeout_check) > 0 && timeout_check[1] != "" && !is.na(timeout_check[1])) {
    .TIMEOUT_CMD_CACHE <<- "timeout"
    return(.TIMEOUT_CMD_CACHE)
  }
  gtimeout_check <- tryCatch(system2("which", "gtimeout", stdout = TRUE, stderr = TRUE), error = function(e) "")
  if (length(gtimeout_check) > 0 && gtimeout_check[1] != "" && !is.na(gtimeout_check[1])) {
    .TIMEOUT_CMD_CACHE <<- "gtimeout"
    return(.TIMEOUT_CMD_CACHE)
  }
  .TIMEOUT_CMD_CACHE <<- ""
  .TIMEOUT_CMD_CACHE
}

build_markov_bgfile <- function(fasta_filepath, bgfile_path, order = 1L) {
  bin <- get_fasta_get_markov_bin()
  if (!nzchar(bin)) {
    warning("fasta-get-markov not found; falling back to motif-file background in FIMO.")
    return(NULL)
  }
  dir.create(dirname(bgfile_path), recursive = TRUE, showWarnings = FALSE)
  cmd <- sprintf("%s -m %d %s %s",
                 shQuote(bin),
                 as.integer(order),
                 shQuote(fasta_filepath),
                 shQuote(bgfile_path))
  status <- system(paste0(cmd, " > /dev/null 2>&1"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (status != 0 || !file.exists(bgfile_path)) {
    warning("Failed to build background model with fasta-get-markov; using motif-file background.")
    return(NULL)
  }
  bgfile_path
}

#' Run FIMO pipeline for motif matching (optimized version with pre-generated FASTA)
#' @param tf_folder Folder containing MEME motif files
#' @param tf_name Transcription factor name
#' @param fimo_folder Output folder for FIMO results
#' @param coords Peak coordinates in format "chr-start-end"
#' @param tg_name Target gene name (for logging)
#' @param p_value_threshold FIMO peak-level binding threshold on site p-values
#' @param fasta_filepath Path to pre-generated FASTA file
#' @return Binary binding vector (TRUE/FALSE for each peak)
run_fimo_pipeline_optimized <- function(tf_folder, tf_name, fimo_folder, coords, tg_name,
                                        p_value_threshold = 0.05,
                                        fasta_filepath,
                                        scan_p_value_threshold = 0.1,
                                        bgfile_path = NULL) {

  # Replace "::" with "__" in TF name for file system usage
  safe_tf_name <- gsub("::", "__", tf_name, fixed = TRUE)

  # Step 1: Prepare FIMO command
  meme_filepath <- file.path(tf_folder, paste0(safe_tf_name, ".meme"))
  fimo_output_dir <- file.path(fimo_folder, paste0("fimo_out_", safe_tf_name))
  fimo_results_filepath <- file.path(fimo_output_dir, "fimo.tsv")

  # Check if MEME file exists
  if (!file.exists(meme_filepath)) {
    return(NULL)
  }

  # Check if FASTA file exists and is not empty
  if (!file.exists(fasta_filepath) || file.info(fasta_filepath)$size == 0) {
    return(NULL)
  }

  # Step 2: Run FIMO with simplified timeout (30 seconds per TF - should be fast)
  timeout_seconds <- 30
  bg_arg <- ""
  if (!is.null(bgfile_path) && file.exists(bgfile_path)) {
    bg_arg <- paste("--bgfile", shQuote(bgfile_path))
  }
  fimo_cmd <- sprintf("%s --thresh %s %s --oc %s %s %s",
                      shQuote(get_fimo_bin()),
                      format(scan_p_value_threshold, scientific = FALSE, trim = TRUE),
                      bg_arg,
                      shQuote(fimo_output_dir),
                      shQuote(meme_filepath),
                      shQuote(fasta_filepath))

  timeout_cmd <- get_timeout_command()
  timeout_available <- nzchar(timeout_cmd)

  # Run FIMO with or without timeout
  if (timeout_available) {
    full_cmd <- sprintf("%s %d %s > /dev/null 2>&1", timeout_cmd, timeout_seconds, fimo_cmd)
    exit_status <- system(full_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Check if timeout occurred (exit status 124 for timeout command)
    if (exit_status == 124) {
      return(NULL)
    } else if (exit_status != 0) {
      return(NULL)
    }
  } else {
    # No timeout command - run directly (should be fast for individual TFs)
    exit_status <- system(paste0(fimo_cmd, " > /dev/null 2>&1"), ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (exit_status != 0) {
      return(NULL)
    }
  }

  # Step 3: Check and read FIMO results
  if (!file.exists(fimo_results_filepath)) {
    return(NULL)
  }

  # Safely read table with error catching
  fimo_results <- tryCatch(
    {read.table(fimo_results_filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)},
    error = function(msg){
      return(NULL)
    })

  if (is.null(fimo_results)){
    return(NULL)
  } else {
    if (nrow(fimo_results) == 0) {
      return(rep(FALSE, length(coords)))
    }

    seq_col <- if ("sequence_name" %in% names(fimo_results)) "sequence_name" else if ("sequence.name" %in% names(fimo_results)) "sequence.name" else NULL
    pval_col <- if ("p.value" %in% names(fimo_results)) "p.value" else if ("p-value" %in% names(fimo_results)) "p-value" else NULL
    if (is.null(seq_col) || is.null(pval_col)) {
      return(rep(FALSE, length(coords)))
    }

    seq_names <- as.character(fimo_results[[seq_col]])
    pvals <- suppressWarnings(as.numeric(fimo_results[[pval_col]]))
    valid <- is.finite(pvals) & !is.na(seq_names) & nzchar(seq_names)
    if (!any(valid)) {
      return(rep(FALSE, length(coords)))
    }

    seq_names <- seq_names[valid]
    pvals <- pvals[valid]
    in_coords <- seq_names %in% coords
    if (!any(in_coords)) {
      return(rep(FALSE, length(coords)))
    }
    seq_names <- seq_names[in_coords]
    pvals <- pvals[in_coords]

    # Peak-level binding call from original FIMO p-values:
    # call bound if any motif-site p-value in that peak is <= threshold.
    binding_vector <- rep(FALSE, length(coords))
    names(binding_vector) <- coords
    pvals_by_peak <- split(pvals, seq_names)
    for (peak_name in names(pvals_by_peak)) {
      pv <- pvals_by_peak[[peak_name]]
      pv <- pv[is.finite(pv)]
      if (length(pv) == 0) next
      binding_vector[peak_name] <- any(pv <= p_value_threshold, na.rm = TRUE)
    }
    binding_vector[is.na(binding_vector)] <- FALSE
    as.logical(binding_vector)
  }
}

#' Run FIMO pipeline for motif matching (original version - kept for compatibility)
#' @param tf_folder Folder containing MEME motif files
#' @param tf_name Transcription factor name
#' @param fimo_folder Output folder for FIMO results
#' @param coords Peak coordinates in format "chr-start-end"
#' @param tg_name Target gene name (for FASTA file naming)
#' @param p_value_threshold FIMO peak-level binding threshold on site p-values
#' @param fasta_output_dir Directory to save FASTA file (default: current directory)
#' @return Binary binding vector (TRUE/FALSE for each peak)
run_fimo_pipeline <- function(tf_folder, tf_name, fimo_folder, coords, tg_name,
                              p_value_threshold = 0.05,
                              fasta_output_dir = ".",
                              scan_p_value_threshold = 0.1,
                              bgfile_path = NULL) {

  # Replace "::" with "__" in TF name for file system usage
  safe_tf_name <- gsub("::", "__", tf_name, fixed = TRUE)

  # Step 1: Prepare the Genomic Coordinates
  gr <- do.call(c, lapply(coords, function(x) {
    parts <- unlist(strsplit(x, "[-]"))
    chr <- parts[1]
    range_start <- as.numeric(parts[2])
    range_end <- as.numeric(parts[3])
    GRanges(seqnames = chr,
            ranges = IRanges(range_start, range_end))
  }))

  # Step 2: Extract DNA sequences
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  names(seqs) <- coords

  # Step 3: Write sequences to FASTA in specified directory
  if (!dir.exists(fasta_output_dir)) {
    dir.create(fasta_output_dir, recursive = TRUE)
  }
  fasta_filepath <- file.path(fasta_output_dir, sprintf("%s_cre_sequences.fa", tg_name))
  writeXStringSet(seqs, filepath = fasta_filepath)

  # Step 4: Run FIMO with timeout (suppress output)
  meme_filepath <- file.path(tf_folder, paste0(safe_tf_name, ".meme"))
  fimo_output_dir <- file.path(fimo_folder, paste0("fimo_out_", safe_tf_name))
  fimo_results_filepath <- file.path(fimo_output_dir, "fimo.tsv")

  # Check if MEME file exists
  if (!file.exists(meme_filepath)) {
    cat("    Warning: MEME file not found:", meme_filepath, "\n")
    return(NULL)
  }

  # Check if FASTA file exists and is not empty
  if (!file.exists(fasta_filepath) || file.info(fasta_filepath)$size == 0) {
    cat("    Warning: FASTA file not found or empty:", fasta_filepath, "\n")
    return(NULL)
  }

  # Check FASTA file size - if too large, it might cause FIMO to hang
  fasta_size_mb <- file.info(fasta_filepath)$size / (1024 * 1024)
  if (fasta_size_mb > 100) {
    cat("    Warning: FASTA file is very large (", round(fasta_size_mb, 1), "MB). FIMO may take a long time or hang.\n", sep = "")
  }

  # Try to use timeout command (available on Linux, macOS with GNU coreutils, or via gtimeout)
  # Set timeout to 5 minutes per TF
  timeout_seconds <- 300
  bg_arg <- ""
  if (!is.null(bgfile_path) && file.exists(bgfile_path)) {
    bg_arg <- paste("--bgfile", shQuote(bgfile_path))
  }
  fimo_cmd <- sprintf("%s --thresh %s %s --oc %s %s %s",
                      shQuote(get_fimo_bin()),
                      format(scan_p_value_threshold, scientific = FALSE, trim = TRUE),
                      bg_arg,
                      shQuote(fimo_output_dir),
                      shQuote(meme_filepath),
                      shQuote(fasta_filepath))

  timeout_cmd <- get_timeout_command()
  timeout_available <- nzchar(timeout_cmd)

  # Run FIMO with or without timeout
  if (timeout_available) {
    full_cmd <- sprintf("%s %d %s > /dev/null 2>&1", timeout_cmd, timeout_seconds, fimo_cmd)
    exit_status <- system(full_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

    # Check if timeout occurred (exit status 124 for timeout command)
    if (exit_status == 124) {
      cat("    Warning: FIMO timed out after", timeout_seconds, "seconds for TF:", tf_name, "\n")
      return(NULL)
    } else if (exit_status != 0) {
      cat("    Warning: FIMO failed (exit status", exit_status, ") for TF:", tf_name, "\n")
      return(NULL)
    }
  } else {
    # No timeout command available - use background process with monitoring
    # Create a wrapper script that runs FIMO
    temp_script <- tempfile(pattern = paste0("fimo_", safe_tf_name, "_"), fileext = ".sh")
    script_content <- c(
      "#!/bin/bash",
      paste0("cd ", shQuote(getwd())),
      fimo_cmd
    )
    writeLines(script_content, temp_script)
    system(paste("chmod +x", shQuote(temp_script)))

    # Run in background and capture PID
    start_time_fimo <- Sys.time()
    bg_output <- system(paste0(temp_script, " > /dev/null 2>&1 & echo $!"), intern = TRUE)
    bg_pid <- as.integer(bg_output[length(bg_output)])

    if (is.na(bg_pid) || bg_pid <= 0) {
      cat("    Warning: Failed to start FIMO process for TF:", tf_name, "\n")
      unlink(temp_script)
      return(NULL)
    }

    # Monitor process with timeout
    process_finished <- FALSE
    waited <- 0
    check_interval <- 5  # Check every 5 seconds

    while (waited < timeout_seconds && !process_finished) {
      Sys.sleep(check_interval)
      waited <- waited + check_interval

      # Check if process is still running using ps
      ps_check <- tryCatch({
        ps_output <- system(paste("ps -p", bg_pid, "2>&1"), intern = TRUE, ignore.stderr = TRUE)
        length(ps_output) > 1 && !any(grepl("No such process", ps_output))
      }, error = function(e) FALSE)

      if (!ps_check) {
        # Process finished
        process_finished <- TRUE
      } else {
        # Still running - check if we've exceeded timeout
        elapsed_fimo <- as.numeric(Sys.time() - start_time_fimo, units = "secs")
        if (elapsed_fimo >= timeout_seconds) {
          # Kill the process
          system(paste("kill -9", bg_pid, "2>/dev/null"), ignore.stdout = TRUE, ignore.stderr = TRUE)
          unlink(temp_script)
          cat("    Warning: FIMO timed out after", timeout_seconds, "seconds for TF:", tf_name, "\n")
          return(NULL)
        }
      }
    }

    # Final check - if still running after timeout, kill it
    if (!process_finished) {
      system(paste("kill -9", bg_pid, "2>/dev/null"), ignore.stdout = TRUE, ignore.stderr = TRUE)
      unlink(temp_script)
      cat("    Warning: FIMO timed out after", timeout_seconds, "seconds for TF:", tf_name, "\n")
      return(NULL)
    }

    # Clean up script
    unlink(temp_script)

    # Wait a moment for file system to sync
    Sys.sleep(1)

    # Check if output file was created
    if (!file.exists(fimo_results_filepath)) {
      cat("    Warning: FIMO output file not created for TF:", tf_name, "\n")
      return(NULL)
    }
  }

  # Step 5: Check and read FIMO results
  # (fimo_results_filepath already defined above)

  # Safely read table with error catching
  fimo_results <- tryCatch(
    {read.table(fimo_results_filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)},
    error = function(msg){
      return(NULL)
    })

  if (is.null(fimo_results)){
    return(NULL)
  } else {
    if (nrow(fimo_results) == 0) {
      return(rep(FALSE, length(coords)))
    }

    seq_col <- if ("sequence_name" %in% names(fimo_results)) "sequence_name" else if ("sequence.name" %in% names(fimo_results)) "sequence.name" else NULL
    pval_col <- if ("p.value" %in% names(fimo_results)) "p.value" else if ("p-value" %in% names(fimo_results)) "p-value" else NULL
    if (is.null(seq_col) || is.null(pval_col)) {
      return(rep(FALSE, length(coords)))
    }

    seq_names <- as.character(fimo_results[[seq_col]])
    pvals <- suppressWarnings(as.numeric(fimo_results[[pval_col]]))
    valid <- is.finite(pvals) & !is.na(seq_names) & nzchar(seq_names)
    if (!any(valid)) {
      return(rep(FALSE, length(coords)))
    }

    seq_names <- seq_names[valid]
    pvals <- pvals[valid]
    in_coords <- seq_names %in% coords
    if (!any(in_coords)) {
      return(rep(FALSE, length(coords)))
    }
    seq_names <- seq_names[in_coords]
    pvals <- pvals[in_coords]

    # Peak-level binding call from original FIMO p-values:
    # call bound if any motif-site p-value in that peak is <= threshold.
    binding_vector <- rep(FALSE, length(coords))
    names(binding_vector) <- coords
    pvals_by_peak <- split(pvals, seq_names)
    for (peak_name in names(pvals_by_peak)) {
      pv <- pvals_by_peak[[peak_name]]
      pv <- pv[is.finite(pv)]
      if (length(pv) == 0) next
      binding_vector[peak_name] <- any(pv <= p_value_threshold, na.rm = TRUE)
    }

    binding_vector[is.na(binding_vector)] <- FALSE
    as.logical(binding_vector)
  }
}

#' Extract TF name from MEME file
#' @param meme_file Path to MEME file
#' @return TF name or NULL
extract_tf_name_from_meme <- function(meme_file) {
  lines <- readLines(meme_file)
  motif_line <- grep("^MOTIF ", lines, value = TRUE)

  if (length(motif_line) == 0) {
    return(NULL)
  }

  tf_name <- strsplit(motif_line, "\\s+")[[1]][3]
  return(tf_name)
}

#' Extract all TF names from JASPAR directory
#' @param jaspar_dir Directory containing MEME files
#' @return Vector of TF names
extract_all_tf_names <- function(jaspar_dir) {
  meme_files <- list.files(jaspar_dir, pattern = "\\.meme$", full.names = TRUE)

  tf_names <- character()
  for (meme_file in meme_files) {
    tf_name <- extract_tf_name_from_meme(meme_file)
    if (!is.null(tf_name)) {
      tf_names <- c(tf_names, tf_name)
    }
  }

  return(unique(tf_names))
}

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths are supplied by the public entry point when this module is sourced.
if (!exists("BASE_DIR", inherits = FALSE)) BASE_DIR <- getwd()
if (!exists("DATA_DIR", inherits = FALSE)) DATA_DIR <- file.path(BASE_DIR, "../data")
if (!exists("ATAC_DIR", inherits = FALSE)) ATAC_DIR <- file.path(DATA_DIR, "atac")
if (!exists("JASPAR_DIR", inherits = FALSE)) JASPAR_DIR <- file.path(DATA_DIR, "JASPAR2024_CORE_pfms_meme")

# Main ArchRProject file (contains all batches)
if (!exists("MAIN_ARCHR_PROJ_FILE", inherits = FALSE)) MAIN_ARCHR_PROJ_FILE <- file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-Save-ArchR-Project.rds")
if (!exists("ATAC_META_FILE", inherits = FALSE)) ATAC_META_FILE <- file.path(ATAC_DIR, "share-grn_results-ArchR_all_with_peaks-ArchR_all_with_peaks-ATAC_P3-P15_meta_D240727.tsv")

# Metacell directory (for finding batches in period)
if (!exists("METACELL_DIR", inherits = FALSE)) METACELL_DIR <- file.path(BASE_DIR, "../explore/explore_results/metacells")

# FIMO files directory (for FASTA sequences)
if (!exists("FIMO_FILES_DIR", inherits = FALSE)) FIMO_FILES_DIR <- file.path(BASE_DIR, "period_motif_matching_fimo_files")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Get all batches belonging to a specific period
#' @param period Period number (e.g., 3, 4, 5)
#' @return Vector of batch names
get_valid_atac_batches_for_period <- function(period) {
  if (!file.exists(ATAC_META_FILE)) {
    stop("ERROR: ATAC metadata file not found: ", ATAC_META_FILE)
  }

  atac_meta <- tryCatch({
    read.csv(ATAC_META_FILE, check.names = FALSE)
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
  # Pattern: mapping_p{period}_{batch}.csv
  pattern <- paste0("^mapping_p", period, "_.*\\.csv$")
  mapping_files <- list.files(METACELL_DIR, pattern = pattern, full.names = FALSE)

  if (length(mapping_files) == 0) {
    return(character(0))
  }

  # Extract batch names from filenames
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

#' Check if period has any batches
#' @param period Period number
#' @return List with available status and message
check_period_availability <- function(period) {
  batches <- get_batches_for_period(period)

  if (length(batches) == 0) {
    return(list(
      available = FALSE,
      message = paste("No batches found for period", period, "in", METACELL_DIR)
    ))
  }

  # Check if main ArchRProject exists
  if (!file.exists(MAIN_ARCHR_PROJ_FILE)) {
    return(list(
      available = FALSE,
      message = paste("Main ArchRProject file not found at:", MAIN_ARCHR_PROJ_FILE)
    ))
  }

  # Verify at least some batches exist in ArchRProject
  tryCatch({
    proj <- readRDS(MAIN_ARCHR_PROJ_FILE)
    available_samples <- unique(proj$Sample)

    # Check if any batches from period exist in project
    matching_batches <- intersect(batches, available_samples)

    if (length(matching_batches) == 0) {
      return(list(
        available = FALSE,
        message = paste("None of the batches for period", period, "found in ArchRProject.",
                       "Period batches:", paste(batches, collapse = ", "))
      ))
    }

    return(list(
      available = TRUE,
      message = paste("Period", period, "has", length(batches), "batches:",
                     paste(batches, collapse = ", "),
                     "(", length(matching_batches), "found in ArchRProject)")
    ))
  }, error = function(e) {
    return(list(
      available = FALSE,
      message = paste("Error loading ArchRProject:", e$message)
    ))
  })
}

# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

#' Run motif matching for target gene at period level
#' @param target_gene Target gene symbol
#' @param region_size Size of region around TSS
#' @param p_value_threshold FIMO p-value threshold
#' @param period Period number
#' @param output_dir Output directory
#' @return List with motif matrix, TFs, peaks, and metadata
run_motif_matching_for_period <- function(target_gene, region_size, p_value_threshold,
                                          period, output_dir) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PERIOD-SPECIFIC MOTIF MATCHING PIPELINE\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Target Gene:", target_gene, "\n")
  cat("Region Size:", region_size, "\n")
  cat("P-value Threshold:", p_value_threshold, "\n")
  cat("Period:", period, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  # Create output directory
  period_output_dir <- file.path(output_dir, paste0("period_", period), target_gene)
  if (!dir.exists(period_output_dir)) {
    dir.create(period_output_dir, recursive = TRUE)
  }

  # Step 1: Check period availability
  cat("Step 1: Checking period availability...\n")
  period_check <- check_period_availability(period)

  if (!period_check$available) {
    error_msg <- paste("ERROR: Period", period, "is not available.\n", period_check$message)
    cat(error_msg, "\n")

    # Write error to file
    error_file <- file.path(output_dir, paste0("ERROR_period_", period, "_", target_gene, ".txt"))
    writeLines(error_msg, error_file)
    cat("Error information saved to:", error_file, "\n")

    stop(error_msg)
  }

  cat("✓", period_check$message, "\n")

  # Get batches in period
  batches <- get_batches_for_period(period)
  cat("  Batches in period", period, ":", paste(batches, collapse = ", "), "\n")

  # Step 2: Load main ArchRProject (peak set is shared across batches)
  cat("\nStep 2: Loading main ArchRProject...\n")

  if (!file.exists(MAIN_ARCHR_PROJ_FILE)) {
    stop("ERROR: ArchRProject file not found at: ", MAIN_ARCHR_PROJ_FILE)
  }

  proj <- readRDS(MAIN_ARCHR_PROJ_FILE)
  cat("✓ Loaded main ArchRProject with", length(proj$cellNames), "cells\n")

  # Get peak set
  cat("\nStep 3: Extracting peak set...\n")
  peak_set <- getPeakSet(proj)

  if (is.null(peak_set) || length(peak_set) == 0) {
    stop("ERROR: Could not load peak set from ArchRProject.")
  }

  cat("✓ Loaded", length(peak_set), "peaks\n")

  # Get target gene coordinates
  cat("\nStep 4: Getting target gene coordinates...\n")
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  gene_info <- getBM(attributes = c("hgnc_symbol", "chromosome_name",
                                   "start_position", "end_position", "strand"),
                     filters = "hgnc_symbol", values = target_gene, mart = mart)

  if (nrow(gene_info) == 0) {
    # Provide helpful error message with common alternative symbols
    alternative_hints <- ""
    if (target_gene == "SUV420H1") {
      alternative_hints <- "\n  Note: SUV420H1 may be stored as 'KMT5B' in some databases."
    } else if (target_gene == "SUV420H2") {
      alternative_hints <- "\n  Note: SUV420H2 may be stored as 'KMT5C' in some databases."
    }

    error_msg <- paste0(
      "ERROR: Target gene '", target_gene, "' not found in Ensembl database.",
      alternative_hints,
      "\n  Please check the gene symbol and use the official HGNC symbol.",
      "\n  The pipeline will stop to prevent hanging."
    )
    cat(error_msg, "\n")

    # Write error to file for better detection by shell scripts
    error_file <- file.path(output_dir, paste0("ERROR_gene_not_found_", target_gene, ".txt"))
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    writeLines(error_msg, error_file)
    cat("Error information saved to:", error_file, "\n")

    stop(error_msg)
  }

  # Handle multiple gene entries (take first)
  if (nrow(gene_info) > 1) {
    cat("  Warning: Multiple entries found for", target_gene, ". Using first entry.\n")
    gene_info <- gene_info[1, ]
  }

  # Calculate TSS and region
  tss <- if (gene_info$strand == 1) gene_info$start_position else gene_info$end_position
  region <- GRanges(seqnames = paste0("chr", gene_info$chromosome_name),
                    ranges = IRanges::IRanges(start = tss - region_size, end = tss + region_size))

  cat("✓ Target gene region:", as.character(region), "\n")
  cat("  Chromosome:", gene_info$chromosome_name, "\n")
  cat("  TSS:", tss, "\n")
  cat("  Strand:", ifelse(gene_info$strand == 1, "+", "-"), "\n")

  # Find nearby peaks
  cat("\nStep 5: Finding nearby peaks...\n")
  overlaps <- findOverlaps(peak_set, region)
  nearby_peaks <- peak_set[queryHits(overlaps)]
  nearby_peaks_coords <- paste(seqnames(nearby_peaks), start(nearby_peaks), end(nearby_peaks), sep = "-")

  cat("✓ Found", length(nearby_peaks), "nearby peaks\n")

  if (length(nearby_peaks) == 0) {
    warning_msg <- paste("WARNING: No peaks found in region for", target_gene, "in period", period)
    cat(warning_msg, "\n")

    # Write warning to file
    warning_file <- file.path(period_output_dir, paste0("WARNING_period_", period, "_", target_gene, ".txt"))
    writeLines(warning_msg, warning_file)

    return(list(
      motif_matrix = NULL,
      tfs = character(0),
      peaks = character(0),
      region = region,
      period = period,
      batches = batches,
      warning = warning_msg
    ))
  }

  # Run motif analysis
  cat("\nStep 6: Running motif matching...\n")

  # Check JASPAR directory
  if (!dir.exists(JASPAR_DIR)) {
    stop("ERROR: JASPAR directory not found at: ", JASPAR_DIR)
  }

  TF_collection <- extract_all_tf_names(JASPAR_DIR)
  cat("✓ Found", length(TF_collection), "TFs in JASPAR database\n")

  # Create FIMO results directory
  fimo_results_dir <- file.path(period_output_dir, "fimo_results")
  if (!dir.exists(fimo_results_dir)) {
    dir.create(fimo_results_dir, recursive = TRUE)
  }

  # Create FIMO files directory for FASTA sequences
  if (!dir.exists(FIMO_FILES_DIR)) {
    dir.create(FIMO_FILES_DIR, recursive = TRUE)
    cat("✓ Created FIMO files directory:", FIMO_FILES_DIR, "\n")
  }

  # Generate FASTA file ONCE before the loop (optimization)
  cat("  Generating FASTA file for", length(nearby_peaks_coords), "peaks...\n")
  gr <- do.call(c, lapply(nearby_peaks_coords, function(x) {
    parts <- unlist(strsplit(x, "[-]"))
    chr <- parts[1]
    range_start <- as.numeric(parts[2])
    range_end <- as.numeric(parts[3])
    GRanges(seqnames = chr,
            ranges = IRanges(range_start, range_end))
  }))
  seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, gr)
  names(seqs) <- nearby_peaks_coords
  fasta_filepath <- file.path(FIMO_FILES_DIR, sprintf("%s_cre_sequences.fa", target_gene))
  writeXStringSet(seqs, filepath = fasta_filepath)
  cat("✓ FASTA file generated:", fasta_filepath, "\n")
  bgfile_path <- file.path(FIMO_FILES_DIR, sprintf("%s_cre_sequences.bg", target_gene))
  bgfile_built <- build_markov_bgfile(
    fasta_filepath = fasta_filepath,
    bgfile_path = bgfile_path,
    order = 1L
  )
  if (!is.null(bgfile_built)) {
    cat("✓ Background model generated with fasta-get-markov:", bgfile_built, "\n")
  } else {
    cat("  Warning: using motif-file background (no external bg model available)\n")
  }

  # Check FASTA file size - if too large, it might cause FIMO to hang
  fasta_size_mb <- file.info(fasta_filepath)$size / (1024 * 1024)
  if (fasta_size_mb > 100) {
    cat("  Warning: FASTA file is large (", round(fasta_size_mb, 1), "MB). FIMO may take longer.\n", sep = "")
  }

  # Run FIMO for each TF
  cat("\n  Running FIMO for", length(TF_collection), "TFs...\n")
  cat("  This may take a while (processing", length(TF_collection), "TFs)...\n")
  motif_library <- list()
  pb <- txtProgressBar(min = 0, max = length(TF_collection), style = 3, char = "=")

  start_time <- Sys.time()

  for (i in seq_along(TF_collection)) {
    tf <- TF_collection[i]

    # Log progress every 50 TFs or for first/last, or every 10 TFs for better visibility
    if (i %% 10 == 0 || i == 1 || i == length(TF_collection)) {
      elapsed <- as.numeric(Sys.time() - start_time, units = "secs")
      if (i > 1 && elapsed > 0) {
        rate <- i / elapsed
        remaining <- (length(TF_collection) - i) / rate
        cat("\n  Progress: TF", i, "of", length(TF_collection),
            "(", round(i/length(TF_collection)*100, 1), "%)",
            "- Elapsed:", round(elapsed/60, 1), "min",
            "- Estimated remaining:", round(remaining/60, 1), "min",
            "- Current TF:", tf, "\n")
      } else {
        cat("\n  Progress: TF", i, "of", length(TF_collection),
            "- Current TF:", tf, "\n")
      }
      flush.console()
    }

    # Run FIMO with timeout protection (using pre-generated FASTA file)
    binding_vector <- tryCatch({
      run_fimo_pipeline_optimized(JASPAR_DIR, tf, fimo_results_dir,
                                  nearby_peaks_coords, target_gene, p_value_threshold,
                                  fasta_filepath = fasta_filepath,
                                  bgfile_path = bgfile_built)
    }, error = function(e) {
      cat("    Error running FIMO for TF", tf, ":", e$message, "\n")
      return(NULL)
    })

    if (!is.null(binding_vector)) {
      motif_library[[tf]] <- binding_vector
    } else {
      # Create zero vector if FIMO failed
      motif_library[[tf]] <- rep(FALSE, length(nearby_peaks_coords))
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  total_time <- as.numeric(Sys.time() - start_time, units = "secs")
  cat("\n  ✓ Completed FIMO for all", length(TF_collection), "TFs in",
      round(total_time/60, 1), "minutes\n")

  # Create motif matrix
  motif_matrix <- do.call(rbind, motif_library)
  colnames(motif_matrix) <- nearby_peaks_coords

  # Convert to binary (0/1) matrix
  motif_matrix <- as.matrix(motif_matrix)
  storage.mode(motif_matrix) <- "integer"

  cat("\n✓ Created motif matrix:", dim(motif_matrix), "\n")
  cat("  TFs with at least one binding:", sum(rowSums(motif_matrix) > 0), "\n")
  cat("  Peaks with at least one TF binding:", sum(colSums(motif_matrix) > 0), "\n")

  # Save results
  cat("\nStep 7: Saving results...\n")

  # Save motif matrix
  motif_matrix_file <- file.path(period_output_dir, paste0(target_gene, "_motif_matrix.csv"))
  write.csv(motif_matrix, motif_matrix_file)
  cat("✓ Saved motif matrix to:", motif_matrix_file, "\n")

  # Save metadata
  metadata <- list(
    target_gene = target_gene,
    period = period,
    batches = batches,
    region_size = region_size,
    p_value_threshold = p_value_threshold,
    tss_coordinates = as.character(region),
    n_peaks = length(nearby_peaks),
    n_tfs = nrow(motif_matrix),
    tf_names = rownames(motif_matrix),
    peak_coordinates = nearby_peaks_coords,
    n_tfs_with_binding = sum(rowSums(motif_matrix) > 0),
    n_peaks_with_binding = sum(colSums(motif_matrix) > 0)
  )

  metadata_file <- file.path(period_output_dir, paste0(target_gene, "_metadata.rds"))
  saveRDS(metadata, metadata_file)
  cat("✓ Saved metadata to:", metadata_file, "\n")

  # Save summary text file
  summary_file <- file.path(period_output_dir, paste0(target_gene, "_summary.txt"))
  summary_lines <- c(
    paste("Period-Specific Motif Matching Pipeline Summary"),
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Period:", period),
    paste("Batches in Period:", paste(batches, collapse = ", ")),
    paste("Region Size:", region_size),
    paste("P-value Threshold:", p_value_threshold),
    paste(""),
    paste("Results:"),
    paste("  Total TFs tested:", nrow(motif_matrix)),
    paste("  TFs with binding:", sum(rowSums(motif_matrix) > 0)),
    paste("  Total peaks in region:", length(nearby_peaks)),
    paste("  Peaks with TF binding:", sum(colSums(motif_matrix) > 0)),
    paste(""),
    paste("Output Directory:", period_output_dir)
  )
  writeLines(summary_lines, summary_file)
  cat("✓ Saved summary to:", summary_file, "\n")

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PERIOD-SPECIFIC MOTIF MATCHING COMPLETED SUCCESSFULLY\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")

  return(list(
    motif_matrix = motif_matrix,
    tfs = rownames(motif_matrix),
    peaks = nearby_peaks_coords,
    region = region,
    period = period,
    batches = batches,
    metadata = metadata
  ))
}
