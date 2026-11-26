# ==============================================================================
# cicero_core_pipeline_fast.R
#
# Speed-oriented additions on top of 03_cicero_core_pipeline.R:
#   1. load_group_assets()      -- reads a group's matrix/coords ONCE
#   2. run_cicero_group_fast()  -- same pipeline as run_cicero_group(), but
#                                   takes preloaded assets (no disk I/O per
#                                   call) and can accept a fixed
#                                   distance_parameter to skip re-estimating
#                                   it on every subsample.
#
# Does not modify 03_cicero_core_pipeline.R -- sources it for the shared
# helpers (peak_chrom, peak_midpoint, classify_pair_type,
# load_genomic_and_peaktype). Source this file instead of / in addition to
# 03_cicero_core_pipeline.R from step3h_fast.R.
# ==============================================================================

source("03_cicero_core_pipeline.R")

# ------------------------------------------------------------------
# load_group_assets(grp_dir)
#
# Reads everything for one group ONCE: the full peak matrix (all cells),
# peak metadata, and the full reduced-coordinate table (UMAP or LSI,
# whichever was exported). Returns a list to be reused across every
# rep/subsample of that group.
# ------------------------------------------------------------------
load_group_assets <- function(grp_dir) {
  grp <- basename(grp_dir)

  required_files <- c("peak_matrix.mtx", "peak_ids.txt", "cell_barcodes.txt", "peak_coords.tsv")
  missing_files <- required_files[!file.exists(file.path(grp_dir, required_files))]
  if (length(missing_files) > 0) stop("[", grp, "] Missing files: ", paste(missing_files, collapse = ", "))

  peak_ids <- readLines(file.path(grp_dir, "peak_ids.txt"))
  peak_ids <- peak_ids[nzchar(peak_ids)]

  cell_ids_all <- readLines(file.path(grp_dir, "cell_barcodes.txt"))
  cell_ids_all <- cell_ids_all[nzchar(cell_ids_all)]

  mat <- Matrix::readMM(file.path(grp_dir, "peak_matrix.mtx"))
  if (inherits(mat, "ngTMatrix")) {
    mat <- Matrix::sparseMatrix(
      i = mat@i + 1L, j = mat@j + 1L, x = rep(1, length(mat@i)),
      dims = dim(mat), repr = "C"
    )
  } else {
    mat <- methods::as(mat, "dgCMatrix")
  }
  if (nrow(mat) != length(peak_ids)) stop("[", grp, "] Matrix rows != peak IDs")
  if (ncol(mat) != length(cell_ids_all)) stop("[", grp, "] Matrix cols != cell IDs")
  rownames(mat) <- peak_ids
  colnames(mat) <- cell_ids_all

  peak_df <- read.delim(file.path(grp_dir, "peak_coords.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
  if (nrow(peak_df) != length(peak_ids)) stop("[", grp, "] Peak metadata rows != peak IDs")
  rownames(peak_df) <- peak_ids
  if (!"gene_short_name" %in% colnames(peak_df)) peak_df$gene_short_name <- rownames(peak_df)

  reduced_type_file <- file.path(grp_dir, "reduced_dim_type.txt")
  reduced_type <- if (file.exists(reduced_type_file)) tolower(trimws(readLines(reduced_type_file)[1])) else "lsi"

  if (reduced_type == "umap" && file.exists(file.path(grp_dir, "umap_coords.tsv"))) {
    coords_full <- read.delim(file.path(grp_dir, "umap_coords.tsv"), row.names = 1, check.names = FALSE)
    if (!all(cell_ids_all %in% rownames(coords_full))) stop("[", grp, "] UMAP coords missing some cells")
  } else {
    lsi_file <- file.path(grp_dir, "lsi_coords.tsv")
    if (!file.exists(lsi_file)) stop("[", grp, "] Missing LSI coordinates")
    coords_full <- read.delim(lsi_file, row.names = 1, check.names = FALSE)
    if (!all(cell_ids_all %in% rownames(coords_full))) stop("[", grp, "] LSI coords missing some cells")
    reduced_type <- "lsi"
  }

  list(grp = grp, mat = mat, peak_df = peak_df, cell_ids_all = cell_ids_all,
       reduced_type = reduced_type, coords_full = coords_full)
}

# ------------------------------------------------------------------
# .build_cicero_cds_from_assets()  (internal helper, used twice below:
# once to get a fixed distance_parameter from the full group, and once
# per-rep on the subsample)
# ------------------------------------------------------------------
.build_cicero_cds_from_assets <- function(assets, cell_ids, k, umap_seed, verbose, vmsg) {
  mat_sub <- assets$mat[, cell_ids, drop = FALSE]
  peak_df <- assets$peak_df

  pd <- methods::new("AnnotatedDataFrame", data = data.frame(row.names = cell_ids))
  fd <- methods::new("AnnotatedDataFrame", data = peak_df)
  cds <- suppressWarnings(monocle::newCellDataSet(
    mat_sub, phenoData = pd, featureData = fd,
    expressionFamily = VGAM::binomialff(), lowerDetectionLimit = 0
  ))

  if (assets$reduced_type == "umap") {
    coords <- assets$coords_full[cell_ids, 1:2, drop = FALSE]
  } else {
    lsi <- assets$coords_full[cell_ids, , drop = FALSE]
    if (any(!is.finite(as.matrix(lsi)))) stop("LSI contains non-finite values")
    if (!requireNamespace("uwot", quietly = TRUE)) stop("uwot is required")
    set.seed(umap_seed)
    n_neighbors <- min(30, nrow(lsi) - 1)
    coords <- as.data.frame(uwot::umap(
      as.matrix(lsi), n_neighbors = n_neighbors, min_dist = 0.3,
      metric = "cosine", n_components = 2, verbose = FALSE
    ))
    rownames(coords) <- cell_ids
  }
  colnames(coords) <- c("UMAP1", "UMAP2")
  if (any(!is.finite(as.matrix(coords)))) stop("UMAP contains non-finite values")

  vmsg("make_cicero_cds()...")
  cicero_cds <- cicero::make_cicero_cds(cds, reduced_coordinates = coords, k = k)
  list(cicero_cds = cicero_cds, mat_sub = mat_sub)
}

# ------------------------------------------------------------------
# estimate_group_distance_parameter()
#
# Runs make_cicero_cds() + estimate_distance_parameter() ONCE using the
# group's FULL cell set (not a subsample). Use its output as
# fixed_distance_parameter in run_cicero_group_fast() to skip
# re-estimating this on every rep. This is an approximation -- it
# assumes the distance parameter doesn't meaningfully change across
# subsamples of the same group -- traded for a large speed gain. Turn
# this off (pass fixed_distance_parameter = NULL to every rep instead)
# if you want the exact original per-rep behavior.
# ------------------------------------------------------------------
estimate_group_distance_parameter <- function(assets, genomic_coords, k = 50,
                                               distance_constraint = 5e5, window = 1e6,
                                               umap_seed = 123, sample_num = 100, verbose = TRUE) {
  vmsg <- function(...) if (verbose) message("[", assets$grp, " | full-group dist param] ", ...)
  built <- .build_cicero_cds_from_assets(assets, assets$cell_ids_all, k, umap_seed, verbose, vmsg)
  vmsg("estimate_distance_parameter() on full group (", ncol(built$mat_sub), " cells)...")
  dpr <- cicero::estimate_distance_parameter(
    built$cicero_cds, window = window, distance_constraint = distance_constraint,
    genomic_coords = genomic_coords, sample_num = sample_num
  )
  dv <- suppressWarnings(as.numeric(unlist(dpr)))
  dv <- dv[is.finite(dv)]
  if (length(dv) == 0) stop("No finite distance parameters for full group ", assets$grp)
  mean(dv)
}

# ------------------------------------------------------------------
# run_cicero_group_fast()
#
# Same behavior/return value as run_cicero_group() in
# 03_cicero_core_pipeline.R, but:
#   - takes `assets` from load_group_assets() instead of grp_dir
#     (no repeated disk reads)
#   - accepts fixed_distance_parameter: if not NULL, skips
#     estimate_distance_parameter() entirely and uses this value
# ------------------------------------------------------------------
run_cicero_group_fast <- function(assets,
                                   cell_subset,
                                   genomic_coords,
                                   peak_type_lookup,
                                   k = 50,
                                   distance_constraint = 5e5,
                                   window = 1e6,
                                   threshold = 0.2,
                                   umap_seed = 123,
                                   fixed_distance_parameter = NULL,
                                   out_prefix = NULL,
                                   verbose = TRUE) {

  grp <- assets$grp
  vmsg <- function(...) if (verbose) message("[", grp, "] ", ...)

  empty_summary <- function() {
    data.frame(group = grp, n_cells = length(cell_subset), n_connections = NA_integer_,
               n_ccans = NA_integer_, mean_ccan_size = NA_real_, median_ccan_size = NA_real_,
               n_peaks_in_ccans = NA_integer_, median_distance_bp = NA_real_, stringsAsFactors = FALSE)
  }

  result <- tryCatch({
    if (!all(cell_subset %in% assets$cell_ids_all)) stop("cell_subset contains barcodes not in this group")

    built <- .build_cicero_cds_from_assets(assets, cell_subset, k, umap_seed, verbose, vmsg)
    cicero_cds <- built$cicero_cds

    if (is.null(fixed_distance_parameter)) {
      vmsg("estimate_distance_parameter()...")
      dpr <- cicero::estimate_distance_parameter(
        cicero_cds, window = window, distance_constraint = distance_constraint,
        genomic_coords = genomic_coords, sample_num = 100
      )
      dv <- suppressWarnings(as.numeric(unlist(dpr)))
      dv <- dv[is.finite(dv)]
      if (length(dv) == 0) stop("No finite distance parameters")
      mean_distance_parameter <- mean(dv)
    } else {
      mean_distance_parameter <- fixed_distance_parameter
      vmsg("Reusing fixed distance_parameter = ", signif(mean_distance_parameter, 6))
    }

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
    chrom1 <- peak_chrom(conns$Peak1); chrom2 <- peak_chrom(conns$Peak2)
    mid1 <- peak_midpoint(conns$Peak1); mid2 <- peak_midpoint(conns$Peak2)
    conns[, distance := ifelse(chrom1 == chrom2, abs(mid1 - mid2), NA_real_)]

    vmsg("generate_ccans()...")
    ccans <- tryCatch(cicero::generate_ccans(conns, coaccess_cutoff_override = threshold), error = function(e) NULL)
    ccans <- if (!is.null(ccans) && nrow(ccans) > 0) data.table::as.data.table(ccans) else data.table::data.table()

    n_ccans <- if (nrow(ccans) > 0) length(unique(ccans$CCAN)) else 0L
    ccan_sizes <- if (nrow(ccans) > 0) ccans[, .N, by = CCAN]$N else numeric(0)

    summary_row <- data.frame(
      group = grp, n_cells = ncol(built$mat_sub), n_connections = nrow(conns), n_ccans = n_ccans,
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
