# ==============================================================================
# 03_cicero_core_pipeline.R
#
# Shared helper: runs the same Cicero pipeline as 02_run_cicero_from_export.R
# (make_cicero_cds -> estimate_distance_parameter -> generate_cicero_models ->
# assemble_connections -> threshold filter -> generate_ccans), but on an
# arbitrary SUBSET of cells from one exported group. This is the one piece
# of logic that step3g (downsample-and-reprocess) and step3h (bootstrap/
# subsample CI) both need, so it lives here once instead of being copied
# into each script.
#
# Not meant to be run directly -- source() it from step3g / step3h.
#
# Environment: same as 02_run_cicero_from_export.R
#   R 4.5.3, monocle 2.x, monocle3 1.4.27, cicero 1.29.2, Matrix 1.7.6,
#   VGAM, stringr, data.table, uwot (only needed for groups exported with
#   reduced_dim_type == "lsi").
# ==============================================================================

suppressPackageStartupMessages({
  library(monocle)
  library(cicero)
  library(Matrix)
  library(VGAM)
  library(stringr)
  library(data.table)
})

# ------------------------------------------------------------------
# Peak helpers (same definitions as step1_cicero_qc.R / step3_cicero_prenatal_postnatal.R)
# ------------------------------------------------------------------
peak_chrom <- function(peak) {
  str_split_fixed(peak, "_", 3)[, 1]
}

peak_midpoint <- function(peak) {
  parts <- str_split_fixed(peak, "_", 3)
  (as.numeric(parts[, 2]) + as.numeric(parts[, 3])) / 2
}

# ------------------------------------------------------------------
# Pair-type classification (verbatim logic from 02_run_cicero_from_export.R)
# ------------------------------------------------------------------
classify_pair_type <- function(peak1_ids, peak2_ids, lookup) {
  t1 <- lookup[peak1_ids]
  t2 <- lookup[peak2_ids]
  t1[is.na(t1) | t1 == ""] <- "Unknown"
  t2[is.na(t2) | t2 == ""] <- "Unknown"

  simplify_type <- function(x) {
    ifelse(x == "Promoter", "Promoter", ifelse(x == "Unknown", "Unknown", "Distal"))
  }

  s1 <- simplify_type(t1)
  s2 <- simplify_type(t2)
  apply(cbind(s1, s2), 1, function(r) paste(sort(r), collapse = "-"))
}

# ------------------------------------------------------------------
# Load the two files shared across every group (once per script run)
# ------------------------------------------------------------------
load_genomic_and_peaktype <- function(export_dir) {
  genomic_coords_file <- file.path(export_dir, "genomic_coords.tsv")
  peak_type_file <- file.path(export_dir, "peak_type_lookup.tsv")
  if (!file.exists(genomic_coords_file)) stop("Missing genomic_coords.tsv: ", genomic_coords_file)
  if (!file.exists(peak_type_file)) stop("Missing peak_type_lookup.tsv: ", peak_type_file)

  genomic_coords <- read.delim(genomic_coords_file, stringsAsFactors = FALSE, check.names = FALSE)

  peak_type_df <- read.delim(peak_type_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("peak_id", "peakType") %in% colnames(peak_type_df))) {
    stop("peak_type_lookup.tsv must contain peak_id and peakType")
  }
  peak_type_lookup <- setNames(peak_type_df$peakType, peak_type_df$peak_id)

  list(genomic_coords = genomic_coords, peak_type_lookup = peak_type_lookup)
}

# ------------------------------------------------------------------
# run_cicero_group()
#
# grp_dir              path to cicero_export/<group>/
# cell_subset           character vector of cell barcodes to use (must be a
#                        subset of that group's cell_barcodes.txt).
#                        NULL = use every cell in the group (full N).
# genomic_coords /
# peak_type_lookup      from load_genomic_and_peaktype()
# k, distance_constraint,
# window, threshold     same Cicero parameters as 02_run_cicero_from_export.R
#                        (defaults below match its CICERO_K / _DISTANCE_
#                        CONSTRAINT / _WINDOW / COACCESS_THRESHOLD)
# umap_seed              seed for uwot::umap() -- only used when the group
#                        was exported as reduced_dim_type == "lsi" and a
#                        fresh 2D embedding has to be computed for the
#                        (sub)set of cells being used here
# out_prefix             if not NULL, writes "<out_prefix>_connections.tsv"
#                        and "<out_prefix>_ccans.tsv" (already filtered to
#                        |coaccess| > threshold), same format as
#                        02_run_cicero_from_export.R's outputs
#
# Returns a list:
#   ok       TRUE/FALSE
#   error    NULL or an error message
#   summary  1-row data.frame: group, n_cells, n_connections, n_ccans,
#            mean_ccan_size, median_ccan_size, n_peaks_in_ccans,
#            median_distance_bp  (matches step1_summary_stats.tsv columns)
#   conns    data.table of filtered connections, or NULL on failure
#   ccans    data.table of CCANs, or NULL on failure
# ------------------------------------------------------------------
run_cicero_group <- function(grp_dir,
                              cell_subset = NULL,
                              genomic_coords,
                              peak_type_lookup,
                              k = 50,
                              distance_constraint = 5e5,
                              window = 1e6,
                              threshold = 0.2,
                              umap_seed = 123,
                              out_prefix = NULL,
                              verbose = TRUE) {

  grp <- basename(grp_dir)
  vmsg <- function(...) if (verbose) message("[", grp, "] ", ...)

  empty_summary <- function() {
    data.frame(
      group = grp,
      n_cells = if (is.null(cell_subset)) NA_integer_ else length(cell_subset),
      n_connections = NA_integer_,
      n_ccans = NA_integer_,
      mean_ccan_size = NA_real_,
      median_ccan_size = NA_real_,
      n_peaks_in_ccans = NA_integer_,
      median_distance_bp = NA_real_,
      stringsAsFactors = FALSE
    )
  }

  result <- tryCatch({

    # ---- required files ----
    required_files <- c("peak_matrix.mtx", "peak_ids.txt", "cell_barcodes.txt", "peak_coords.tsv")
    missing_files <- required_files[!file.exists(file.path(grp_dir, required_files))]
    if (length(missing_files) > 0) stop("Missing files: ", paste(missing_files, collapse = ", "))

    peak_ids <- readLines(file.path(grp_dir, "peak_ids.txt"))
    peak_ids <- peak_ids[nzchar(peak_ids)]

    cell_ids_all <- readLines(file.path(grp_dir, "cell_barcodes.txt"))
    cell_ids_all <- cell_ids_all[nzchar(cell_ids_all)]

    # ---- read + convert peak matrix (same ngTMatrix handling as the main run) ----
    mat <- Matrix::readMM(file.path(grp_dir, "peak_matrix.mtx"))
    if (inherits(mat, "ngTMatrix")) {
      mat <- Matrix::sparseMatrix(
        i = mat@i + 1L, j = mat@j + 1L, x = rep(1, length(mat@i)),
        dims = dim(mat), repr = "C"
      )
    } else {
      mat <- methods::as(mat, "dgCMatrix")
    }

    if (nrow(mat) != length(peak_ids)) stop("Matrix rows (", nrow(mat), ") != peak IDs (", length(peak_ids), ")")
    if (ncol(mat) != length(cell_ids_all)) stop("Matrix cols (", ncol(mat), ") != cell IDs (", length(cell_ids_all), ")")

    rownames(mat) <- peak_ids
    colnames(mat) <- cell_ids_all

    # ---- subset cells: this is the whole point of this function ----
    if (is.null(cell_subset)) {
      cell_ids <- cell_ids_all
    } else {
      if (!all(cell_subset %in% cell_ids_all)) stop("cell_subset contains barcodes not present in this group")
      cell_ids <- cell_subset
    }
    mat <- mat[, cell_ids, drop = FALSE]
    vmsg("Using ", ncol(mat), " / ", length(cell_ids_all), " cells")

    # ---- peak metadata ----
    peak_df <- read.delim(file.path(grp_dir, "peak_coords.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(peak_df) != length(peak_ids)) stop("Peak metadata rows != peak IDs")
    rownames(peak_df) <- peak_ids
    if (!"gene_short_name" %in% colnames(peak_df)) peak_df$gene_short_name <- rownames(peak_df)

    pd <- methods::new("AnnotatedDataFrame", data = data.frame(row.names = cell_ids))
    fd <- methods::new("AnnotatedDataFrame", data = peak_df)

    cds <- suppressWarnings(monocle::newCellDataSet(
      mat, phenoData = pd, featureData = fd,
      expressionFamily = VGAM::binomialff(), lowerDetectionLimit = 0
    ))

    # ---- reduced coordinates, subset to the same cell_ids used above ----
    reduced_type_file <- file.path(grp_dir, "reduced_dim_type.txt")
    reduced_type <- if (file.exists(reduced_type_file)) {
      tolower(trimws(readLines(reduced_type_file)[1]))
    } else {
      "lsi"
    }

    if (reduced_type == "umap" && file.exists(file.path(grp_dir, "umap_coords.tsv"))) {

      # Reuse the exported (global) UMAP, subset down to the cells being used here --
      # this mirrors 01_export_for_cicero.R subsetting the project-wide UMAP per group.
      umap_coords <- read.delim(file.path(grp_dir, "umap_coords.tsv"), row.names = 1, check.names = FALSE)
      if (!all(cell_ids %in% rownames(umap_coords))) stop("UMAP coords missing some cells")
      umap_coords <- umap_coords[cell_ids, 1:2, drop = FALSE]

    } else {

      lsi_file <- file.path(grp_dir, "lsi_coords.tsv")
      if (!file.exists(lsi_file)) stop("Missing LSI coordinates")
      lsi <- read.delim(lsi_file, row.names = 1, check.names = FALSE)
      if (!all(cell_ids %in% rownames(lsi))) stop("LSI coords missing some cells")
      lsi <- lsi[cell_ids, , drop = FALSE]
      if (any(!is.finite(as.matrix(lsi)))) stop("LSI contains non-finite values")

      if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required")

      # Fit UMAP fresh on this (sub)set of cells -- matches what
      # 02_run_cicero_from_export.R does per-group when no ArchR UMAP exists,
      # so a downsampled/bootstrapped rep gets an embedding fit the same
      # way a genuinely smaller group's would be.
      set.seed(umap_seed)
      n_neighbors <- min(30, nrow(lsi) - 1)
      umap_matrix <- uwot::umap(
        as.matrix(lsi), n_neighbors = n_neighbors, min_dist = 0.3,
        metric = "cosine", n_components = 2, verbose = FALSE
      )
      umap_coords <- as.data.frame(umap_matrix)
      rownames(umap_coords) <- cell_ids
      rm(lsi)
    }

    colnames(umap_coords) <- c("UMAP1", "UMAP2")
    if (any(!is.finite(as.matrix(umap_coords)))) stop("UMAP contains non-finite values")

    # ---- Cicero pipeline proper ----
    vmsg("make_cicero_cds()...")
    cicero_cds <- cicero::make_cicero_cds(cds, reduced_coordinates = umap_coords, k = k)

    vmsg("estimate_distance_parameter()...")
    dist_param_result <- cicero::estimate_distance_parameter(
      cicero_cds, window = window, distance_constraint = distance_constraint,
      genomic_coords = genomic_coords, sample_num = 100
    )
    dist_values <- suppressWarnings(as.numeric(unlist(dist_param_result)))
    dist_values <- dist_values[is.finite(dist_values)]
    if (length(dist_values) == 0) stop("No finite distance parameters")
    mean_distance_parameter <- mean(dist_values)

    vmsg("generate_cicero_models()...")
    cicero_models <- cicero::generate_cicero_models(
      cicero_cds, distance_parameter = mean_distance_parameter,
      window = window, genomic_coords = genomic_coords
    )

    vmsg("assemble_connections()...")
    conns <- cicero::assemble_connections(cicero_models, silent = TRUE)
    if (is.null(conns) || nrow(conns) == 0) stop("assemble_connections returned no rows")

    keep <- !is.na(conns$coaccess) & abs(conns$coaccess) > threshold
    conns <- conns[keep, , drop = FALSE]
    if (nrow(conns) == 0) stop("No connections pass threshold ", threshold)

    conns <- data.table::as.data.table(conns)
    conns[, pair_type := classify_pair_type(Peak1, Peak2, peak_type_lookup)]
    chrom1 <- peak_chrom(conns$Peak1)
    chrom2 <- peak_chrom(conns$Peak2)
    mid1 <- peak_midpoint(conns$Peak1)
    mid2 <- peak_midpoint(conns$Peak2)
    conns[, distance := ifelse(chrom1 == chrom2, abs(mid1 - mid2), NA_real_)]

    vmsg("generate_ccans()...")
    ccans <- tryCatch(
      cicero::generate_ccans(conns, coaccess_cutoff_override = threshold),
      error = function(e) NULL
    )
    ccans <- if (!is.null(ccans) && nrow(ccans) > 0) data.table::as.data.table(ccans) else data.table::data.table()

    n_ccans <- if (nrow(ccans) > 0) length(unique(ccans$CCAN)) else 0L
    ccan_sizes <- if (nrow(ccans) > 0) ccans[, .N, by = CCAN]$N else numeric(0)

    summary_row <- data.frame(
      group = grp,
      n_cells = ncol(mat),
      n_connections = nrow(conns),
      n_ccans = n_ccans,
      mean_ccan_size = if (length(ccan_sizes) > 0) round(mean(ccan_sizes), 3) else NA_real_,
      median_ccan_size = if (length(ccan_sizes) > 0) median(ccan_sizes) else NA_real_,
      n_peaks_in_ccans = nrow(ccans),
      median_distance_bp = suppressWarnings(median(conns$distance, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )

    if (!is.null(out_prefix)) {
      dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
      data.table::fwrite(conns, paste0(out_prefix, "_connections.tsv"), sep = "\t")
      data.table::fwrite(ccans, paste0(out_prefix, "_ccans.tsv"), sep = "\t")
    }

    list(ok = TRUE, error = NULL, summary = summary_row, conns = conns, ccans = ccans)

  }, error = function(e) {
    list(ok = FALSE, error = conditionMessage(e), summary = empty_summary(), conns = NULL, ccans = NULL)
  })

  if (isFALSE(result$ok)) vmsg("FAILED: ", result$error)
  result
}
