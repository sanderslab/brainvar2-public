#!/usr/bin/env Rscript

# ============================================================
# 01_run_cellchat_per_sample_h5ad_FINAL.R
#
# Final per-sample CellChat pipeline for the BrainVar h5ad.
# Designed for a small downstream AHI1 comparison.
#
# IMPORTANT FIXES
#   1. Uses only the homozygous-risk list:
#        AHI1_rs13197384_risk_homo_samples.txt
#   2. Correctly parses CellChat ligand complexes such as:
#        Ach-CHAT_SLC5A7  -> CHAT + SLC5A7
#        Glu-SLC17A7_GLS  -> SLC17A7 + GLS
#   3. Keeps Ensembl -> HGNC mapping explicit and saves diagnostics.
#   4. Aligns metadata to expression columns by cell_id exactly.
#   5. Primary communication model = truncatedMean, trim=0.10.
#   6. Optional triMean sensitivity can be enabled without changing
#      the primary analysis.
#   7. Saves raw LR probabilities, pathway probabilities, p-values,
#      detected pathways, target-pathway diagnostics, and CellChat RDS.
# ============================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(future)
  library(Matrix)
  library(rhdf5)
  library(reticulate)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

h5ad_file <- "/mnt/s3/BrainVar_sestan/BrainVar2_allperiods_merged_RNA.h5ad"
meta_file <- "/mnt/s3/BrainVar_sestan/BrainVar2_allperiods_merged_RNA_meta.csv"

# Use the FULL intended cohorts. Do not auto-select a smaller file.
healthy_fp <- "AHI1_rs13197384_healthy_samples.txt"
risk_fp    <- "AHI1_rs13197384_risk_homo_samples.txt"

SAMPLE_COL <- "library"
BROAD_COL <- "majorclass_v1"
SUBTYPE_COL <- "celltype_v1"
EXPRESSION_LAYER <- "lognorm"

MIN_CELLS <- 10L

# Primary model.
PRIMARY_MODEL_NAME <- "truncatedMean10"
PRIMARY_TYPE <- "truncatedMean"
PRIMARY_TRIM <- 0.10

# Methods sensitivity analysis: run triMean in addition to the primary model.
RUN_TRIMEAN_SENSITIVITY <- TRUE

POPULATION_SIZE <- FALSE
RAW_USE <- TRUE
NBOOT <- 100L
SEED <- 1234L

CellChatDB.use <- CellChatDB.human

future::plan("sequential")
set.seed(SEED)

# ------------------------------------------------------------
# OUTPUTS
# ------------------------------------------------------------

base_out_dir <- "cellchat_results_01"
cellchat_out_dir <- file.path(base_out_dir, "cellchat_objects")
results_out_dir <- file.path(base_out_dir, "results")
gene_map_out_dir <- file.path(base_out_dir, "gene_maps")
log_out_dir <- file.path(base_out_dir, "logs")
config_out_dir <- file.path(base_out_dir, "config")

for (d in c(base_out_dir, cellchat_out_dir, results_out_dir,
            gene_map_out_dir, log_out_dir, config_out_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

models <- list(
  truncatedMean10 = list(type = PRIMARY_TYPE, trim = PRIMARY_TRIM)
)
if (RUN_TRIMEAN_SENSITIVITY) {
  models$triMean <- list(type = "triMean", trim = 0)
}

for (model_name in names(models)) {
  for (group_definition in c("broadtype", "subtype4")) {
    dir.create(file.path(cellchat_out_dir, model_name, group_definition), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(results_out_dir, model_name, group_definition), recursive = TRUE, showWarnings = FALSE)
  }
}

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

safe_write_csv <- function(x, file) {
  data.table::fwrite(as.data.frame(x), file = file)
}

array_to_df <- function(arr, feature_name, value_name = "value") {
  if (is.null(arr)) return(data.frame())
  dn <- dimnames(arr)
  if (length(dn) != 3L) stop("Expected 3D array, got: ", paste(dim(arr), collapse = " x "))
  names(dn) <- c("Sender", "Receiver", feature_name)
  dimnames(arr) <- dn
  as.data.frame.table(arr, responseName = value_name, stringsAsFactors = FALSE)
}

# CellChat stores complex ligands in this database with a family prefix.
# Examples: Ach-CHAT_SLC5A7 and Glu-SLC17A7_GLS.
# Remove the family prefix, then split the complex into individual genes.
extract_ligand_genes <- function(x) {
  x <- as.character(x)
  x <- sub("^[^-]+-", "", x)
  parts <- unlist(strsplit(x, "_", fixed = TRUE), use.names = FALSE)
  parts <- parts[!is.na(parts) & nzchar(parts)]
  unique(parts)
}

save_sample_exports <- function(cellchat, sample_id, model_name, group_definition, out_res_dir) {
  # Full LR probabilities and p-values.
  prob_df <- array_to_df(cellchat@net$prob, "LR", "prob")
  pval_df <- array_to_df(cellchat@net$pval, "LR", "pval")

  if (nrow(prob_df) > 0) {
    lr_df <- merge(prob_df, pval_df,
                   by = c("Sender", "Receiver", "LR"), all = TRUE)
    lr_df$Sample <- sample_id
    lr_df$Model <- model_name
    lr_df$GroupDefinition <- group_definition
    safe_write_csv(lr_df, file.path(out_res_dir, paste0(sample_id, "_LR_all.csv")))

    target_lr <- lr_df[
      grepl("Glutamate|Glu-|GRIA|GRIK|GRM|SLC17|SLC1A|GLS|GLS2|Acetylcholine|Ach-|CHAT|SLC5A7|SLC18A3|CHRN|CHRNA|CHRM",
            lr_df$LR, ignore.case = TRUE),
      , drop = FALSE
    ]
    if (nrow(target_lr) > 0) {
      safe_write_csv(target_lr, file.path(out_res_dir, paste0(sample_id, "_DIAG_glutamate_cholinergic_LR.csv")))
    }
  }

  # Full pathway probabilities and p-values.
  if (!is.null(cellchat@netP$prob)) {
    path_prob <- array_to_df(cellchat@netP$prob, "Pathway", "prob")
    path_pval <- if (!is.null(cellchat@netP$pval)) array_to_df(cellchat@netP$pval, "Pathway", "pval") else data.frame()

    if (nrow(path_prob) > 0) {
      if (nrow(path_pval) > 0) {
        path_df <- merge(path_prob, path_pval,
                         by = c("Sender", "Receiver", "Pathway"), all = TRUE)
      } else {
        path_df <- path_prob
        path_df$pval <- NA_real_
      }
      path_df$Sample <- sample_id
      path_df$Model <- model_name
      path_df$GroupDefinition <- group_definition
      safe_write_csv(path_df, file.path(out_res_dir, paste0(sample_id, "_pathway_all.csv")))

      target_path <- path_df[
        grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline",
              path_df$Pathway, ignore.case = TRUE),
        , drop = FALSE
      ]
      if (nrow(target_path) > 0) {
        safe_write_csv(target_path, file.path(out_res_dir, paste0(sample_id, "_DIAG_glutamate_cholinergic_pathway.csv")))
      }
    }
  }

  # Pathways explicitly represented by the object.
  pathways <- if (!is.null(cellchat@netP$pathways)) cellchat@netP$pathways else character(0)
  safe_write_csv(
    data.frame(Sample = sample_id,
               Model = model_name,
               GroupDefinition = group_definition,
               Pathway = pathways,
               stringsAsFactors = FALSE),
    file.path(out_res_dir, paste0(sample_id, "_pathways_detected.csv"))
  )

  # CellChat's own filtered communication tables.
  lr_sig <- tryCatch(subsetCommunication(cellchat), error = function(e) data.frame())
  if (nrow(lr_sig) > 0) safe_write_csv(lr_sig, file.path(out_res_dir, paste0(sample_id, "_LR_subsetCommunication.csv")))

  path_sig <- tryCatch(subsetCommunication(cellchat, slot.name = "netP"), error = function(e) data.frame())
  if (nrow(path_sig) > 0) safe_write_csv(path_sig, file.path(out_res_dir, paste0(sample_id, "_pathway_subsetCommunication.csv")))
}

# ------------------------------------------------------------
# INPUT CHECKS
# ------------------------------------------------------------

for (f in c(h5ad_file, meta_file, healthy_fp, risk_fp)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

healthy_samples <- unique(trimws(readLines(healthy_fp)))
risk_samples <- unique(trimws(readLines(risk_fp)))
healthy_samples <- healthy_samples[nzchar(healthy_samples)]
risk_samples <- risk_samples[nzchar(risk_samples)]

if (length(intersect(healthy_samples, risk_samples)) > 0) {
  stop("Samples appear in BOTH healthy and risk lists: ",
       paste(intersect(healthy_samples, risk_samples), collapse = ", "))
}

cat("Healthy samples:", length(healthy_samples), "\n")
cat("Risk samples:", length(risk_samples), "\n")

# ------------------------------------------------------------
# METADATA / H5AD ALIGNMENT
# ------------------------------------------------------------

needed_meta_cols <- c("cell_id", "library", "majorclass_v1", "celltype_v1")
meta_csv <- data.table::fread(meta_file, select = needed_meta_cols, showProgress = TRUE)
meta_csv <- as.data.frame(meta_csv, stringsAsFactors = FALSE)
meta_csv$cell_id <- as.character(meta_csv$cell_id)
meta_csv$library <- as.character(meta_csv$library)

h5ad_cells <- as.character(rhdf5::h5read(h5ad_file, "/obs/cell_id"))

meta_match <- match(h5ad_cells, meta_csv$cell_id)
valid_match <- !is.na(meta_match)

h5_meta <- data.frame(
  cell_id = h5ad_cells,
  library = NA_character_,
  majorclass_v1 = NA_character_,
  celltype_v1 = NA_character_,
  stringsAsFactors = FALSE
)

h5_meta$library[valid_match] <- meta_csv$library[meta_match[valid_match]]
h5_meta$majorclass_v1[valid_match] <- meta_csv$majorclass_v1[meta_match[valid_match]]
h5_meta$celltype_v1[valid_match] <- meta_csv$celltype_v1[meta_match[valid_match]]

h5_meta$AHI1_status <- NA_character_
h5_meta$AHI1_status[h5_meta$library %in% healthy_samples] <- "healthy"
h5_meta$AHI1_status[h5_meta$library %in% risk_samples] <- "risk"

available_libraries <- unique(h5_meta$library[!is.na(h5_meta$library)])
all_samples <- intersect(unique(c(healthy_samples, risk_samples)), available_libraries)

sample_status <- data.frame(
  Sample = all_samples,
  AHI1_status = ifelse(all_samples %in% risk_samples, "risk", "healthy"),
  stringsAsFactors = FALSE
)
safe_write_csv(sample_status, file.path(config_out_dir, "samples_used_by_h5ad.csv"))

cat("Samples present in h5ad:", length(all_samples), "\n")
cat("Healthy present in h5ad:", sum(sample_status$AHI1_status == "healthy"), "\n")
cat("Risk present in h5ad:", sum(sample_status$AHI1_status == "risk"), "\n")

# ------------------------------------------------------------
# GENE MAPPING
# ------------------------------------------------------------

gene_ids <- as.character(rhdf5::h5read(h5ad_file, "/var/_index"))
gene_ids_clean <- vapply(strsplit(gene_ids, ".", fixed = TRUE), function(x) x[1], character(1))

gene_map_file <- file.path(gene_map_out_dir, "ensembl_to_symbol.rds")

if (file.exists(gene_map_file)) {
  gene_map <- readRDS(gene_map_file)
} else {
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = unique(gene_ids_clean),
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  map_df <- map_df[
    !is.na(map_df$SYMBOL) & nzchar(map_df$SYMBOL),
    , drop = FALSE
  ]
  map_df <- map_df[
    !duplicated(map_df$ENSEMBL),
    , drop = FALSE
  ]
  gene_map <- map_df
  saveRDS(gene_map, gene_map_file)
}

symbol_match <- match(gene_ids_clean, gene_map$ENSEMBL)
gene_symbols_raw <- gene_map$SYMBOL[symbol_match]
missing_symbol <- is.na(gene_symbols_raw) | !nzchar(gene_symbols_raw)
gene_symbols_raw[missing_symbol] <- gene_ids_clean[missing_symbol]

gene_symbols <- make.unique(gene_symbols_raw)

# ------------------------------------------------------------
# BACKED ANNDATA ACCESS
# ------------------------------------------------------------

anndata_py <- reticulate::import("anndata", convert = FALSE)
numpy_py <- reticulate::import("numpy", convert = FALSE)
adata <- anndata_py$read_h5ad(h5ad_file, backed = "r")

reticulate::py_run_string(
'
def extract_sparse_rows(adata, rows_zero_based, layer_name):
    x = adata.layers[layer_name][rows_zero_based, :]
    x = x.transpose().tocsc()
    return {
        "data": x.data,
        "indices": x.indices,
        "indptr": x.indptr,
        "shape": x.shape
    }
'
)

py_main <- reticulate::import("__main__", convert = FALSE)
extract_sparse_rows <- py_main$extract_sparse_rows

read_one_sample_matrix <- function(h5_rows) {
  rows0 <- as.integer(h5_rows - 1L)
  rows_py <- numpy_py$array(rows0, dtype = "int64")
  sp <- extract_sparse_rows(adata, rows_py, EXPRESSION_LAYER)

  sp_data <- reticulate::py_to_r(reticulate::py_get_item(sp, "data"))
  sp_i <- reticulate::py_to_r(reticulate::py_get_item(sp, "indices"))
  sp_p <- reticulate::py_to_r(reticulate::py_get_item(sp, "indptr"))
  sp_shape <- reticulate::py_to_r(reticulate::py_get_item(sp, "shape"))

  expr <- Matrix::sparseMatrix(
    i = as.integer(sp_i) + 1L,
    p = as.integer(sp_p),
    x = as.numeric(sp_data),
    dims = c(sp_shape[[1]], sp_shape[[2]]),
    repr = "C"
  )

  rownames(expr) <- gene_symbols
  colnames(expr) <- h5ad_cells[h5_rows]
  expr
}

# ------------------------------------------------------------
# GROUPING
# ------------------------------------------------------------

make_subtype4 <- function(meta) {
  meta$cci_group <- NA_character_
  meta$cci_group[startsWith(meta$majorclass_v1, "ExN")] <- "ExN"
  meta$cci_group[startsWith(meta$majorclass_v1, "InN")] <- "InN"
  meta$cci_group[meta$majorclass_v1 %in% c("IT-Immature", "NonIT-Immature")] <- "ExN"
  meta$cci_group[meta$celltype_v1 == "Micro"] <- "Microglia"
  meta$cci_group[meta$celltype_v1 == "OPC"] <- "OPC"
  meta
}

make_broadtype5 <- function(meta) {
  meta$cci_group <- NA_character_
  meta$cci_group[startsWith(meta$majorclass_v1, "ExN")] <- "ExN"
  meta$cci_group[startsWith(meta$majorclass_v1, "InN")] <- "InN"
  meta$cci_group[meta$majorclass_v1 %in% c("IT-Immature", "NonIT-Immature")] <- "ExN"
  meta$cci_group[meta$majorclass_v1 == "Glia"] <- "Glia"
  meta$cci_group[meta$majorclass_v1 == "NNC"] <- "NNC"
  meta
}

# ------------------------------------------------------------
# ONE SAMPLE / ONE MODEL
# ------------------------------------------------------------

run_one <- function(sample_id, group_definition, model_name, model_cfg) {
  cat("\n============================================================\n")
  cat("Sample:", sample_id, "| Group:", group_definition, "| Model:", model_name, "\n")
  cat("============================================================\n")

  out_obj_dir <- file.path(cellchat_out_dir, model_name, group_definition)
  out_res_dir <- file.path(results_out_dir, model_name, group_definition)
  outfile <- file.path(
    out_obj_dir,
    paste0("CCI_CellChat_", group_definition, "_", sample_id, ".rds")
  )

  if (file.exists(outfile)) {
    message("Already exists; skipping: ", outfile)
    return(invisible(NULL))
  }

  base_keep <- h5_meta$library == sample_id & !is.na(h5_meta$AHI1_status)

  if (group_definition == "subtype4") {
    sm <- make_subtype4(h5_meta[base_keep, , drop = FALSE])
    group_levels <- c("ExN", "InN", "Microglia", "OPC")
  } else {
    sm <- make_broadtype5(h5_meta[base_keep, , drop = FALSE])
    group_levels <- c("ExN", "InN", "Glia", "NNC")
  }

  sm <- sm[!is.na(sm$cci_group) & nzchar(sm$cci_group), , drop = FALSE]
  if (nrow(sm) == 0) {
    message("No assigned cells for sample.")
    return(invisible(NULL))
  }

  rownames(sm) <- sm$cell_id
  sm$cellchat_group <- factor(sm$cci_group, levels = group_levels)
  sm$samples <- factor(sample_id)

  counts <- table(sm$cellchat_group)
  safe_write_csv(
    data.frame(
      Sample = sample_id,
      Model = model_name,
      GroupDefinition = group_definition,
      group = names(counts),
      n_cells = as.integer(counts),
      stringsAsFactors = FALSE
    ),
    file.path(out_res_dir, paste0(sample_id, "_cell_counts.csv"))
  )

  valid_groups <- names(counts[counts >= MIN_CELLS])
  if (length(valid_groups) < 2) {
    message("Fewer than 2 groups with >= ", MIN_CELLS, " cells; skipping.")
    return(invisible(NULL))
  }

  keep_group <- sm$cellchat_group %in% valid_groups
  sm <- sm[keep_group, , drop = FALSE]
  sm$cellchat_group <- droplevels(sm$cellchat_group)
  rownames(sm) <- sm$cell_id

  h5_rows <- match(rownames(sm), h5ad_cells)
  if (anyNA(h5_rows)) stop("Some selected cells could not be mapped back to h5ad rows.")

  expr <- read_one_sample_matrix(h5_rows)

  # Robust metadata alignment after expression extraction.
  sm <- sm[match(colnames(expr), rownames(sm)), , drop = FALSE]
  rownames(sm) <- sm$cell_id
  if (!identical(colnames(expr), rownames(sm))) {
    stop("Expression / metadata alignment failed for ", sample_id)
  }

  # ----------------------------------------------------------
  # CRITICAL: parse CellChat database complex ligands into genes.
  # ----------------------------------------------------------
  ligand_genes <- extract_ligand_genes(CellChatDB.use$interaction$ligand)
  receptor_genes <- unique(as.character(CellChatDB.use$interaction$receptor))
  receptor_genes <- receptor_genes[!is.na(receptor_genes) & nzchar(receptor_genes)]
  cellchat_genes <- unique(c(ligand_genes, receptor_genes))

  keep_genes <- rownames(expr) %in% cellchat_genes
  expr <- expr[keep_genes, , drop = FALSE]

  # Save sample-specific target gene retention diagnostics.
  target_genes <- c(
    "CHAT", "SLC5A7", "SLC18A3",
    "SLC17A6", "SLC17A7", "SLC17A8", "GLS", "GLS2",
    "GRIA3", "GRIK3", "GRIK4",
    "CHRNA7", "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"
  )

  safe_write_csv(
    data.frame(
      Sample = sample_id,
      Model = model_name,
      GroupDefinition = group_definition,
      Gene = target_genes,
      PresentInExpression = target_genes %in% rownames(expr),
      stringsAsFactors = FALSE
    ),
    file.path(out_res_dir, paste0(sample_id, "_target_gene_retention.csv"))
  )

  if (nrow(expr) < 50) warning("Very few CellChat genes retained for ", sample_id)

  cellchat <- createCellChat(
    object = expr,
    meta = sm,
    group.by = "cellchat_group",
    do.sparse = TRUE
  )
  cellchat@DB <- CellChatDB.use

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat, min.cells = MIN_CELLS)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # Save expression/filter diagnostics for target genes.
  var_features <- if (!is.null(cellchat@var.features$features)) cellchat@var.features$features else character(0)
  safe_write_csv(
    data.frame(
      Sample = sample_id,
      Model = model_name,
      GroupDefinition = group_definition,
      Gene = target_genes,
      InDataSignaling = target_genes %in% rownames(cellchat@data.signaling),
      OverExpressed = target_genes %in% var_features,
      stringsAsFactors = FALSE
    ),
    file.path(out_res_dir, paste0(sample_id, "_target_gene_CellChat_diagnostic.csv"))
  )

  # Save LRsig before communication probability calculation.
  if (!is.null(cellchat@LR$LRsig) && nrow(cellchat@LR$LRsig) > 0) {
    safe_write_csv(
      cellchat@LR$LRsig,
      file.path(out_res_dir, paste0(sample_id, "_LRsig.csv"))
    )
  }

  cellchat <- computeCommunProb(
    cellchat,
    type = model_cfg$type,
    trim = model_cfg$trim,
    raw.use = RAW_USE,
    population.size = POPULATION_SIZE,
    nboot = NBOOT,
    seed.use = SEED
  )

  # Save probability/p-value arrays BEFORE communication filtering.
  raw_prob <- array_to_df(cellchat@net$prob, "LR", "prob")
  raw_pval <- array_to_df(cellchat@net$pval, "LR", "pval")

  if (nrow(raw_prob) > 0) {
    raw_lr <- merge(raw_prob, raw_pval,
                    by = c("Sender", "Receiver", "LR"), all = TRUE)
    raw_lr$Sample <- sample_id
    raw_lr$Model <- model_name
    raw_lr$GroupDefinition <- group_definition
    safe_write_csv(raw_lr, file.path(out_res_dir, paste0(sample_id, "_LR_raw_before_filter.csv")))
  }

  cellchat <- filterCommunication(cellchat, min.cells = MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  save_sample_exports(cellchat, sample_id, model_name, group_definition, out_res_dir)
  saveRDS(cellchat, outfile)

  cat("SAVED:", outfile, "\n")

  rm(expr, sm, cellchat)
  gc()
  invisible(NULL)
}

# ------------------------------------------------------------
# RUN / REPAIR MISSING CELLCHAT OBJECTS
# ------------------------------------------------------------

# This section uses ALL samples present in the healthy/risk lists
# and verified to be present in the h5ad file.
# Existing CellChat objects are preserved and skipped by run_one().
# Missing objects are created automatically.

cat("\n============================================================\n")
cat("CELLCHAT OBJECT REPAIR / GENERATION\n")
cat("============================================================\n")
cat("Healthy samples in list: ", length(healthy_samples), "\n")
cat("Risk samples in list:    ", length(risk_samples), "\n")
cat("Samples present in h5ad: ", length(all_samples), "\n")

# Explicitly report missing objects before starting.
for (model_name in names(models)) {
  for (group_definition in c("subtype4", "broadtype")) {

    obj_dir <- file.path(
      cellchat_out_dir,
      model_name,
      group_definition
    )

    existing_files <- list.files(
      obj_dir,
      pattern = "^CCI_CellChat_.*\\.rds$",
      full.names = FALSE
    )

    existing_ids <- sub(
      paste0("^CCI_CellChat_", group_definition, "_"),
      "",
      existing_files
    )
    existing_ids <- sub("\\.rds$", "", existing_ids)

    missing_ids <- setdiff(all_samples, existing_ids)

    cat("\n------------------------------------------------------------\n")
    cat("Model:", model_name, "| Group:", group_definition, "\n")
    cat("Existing objects:", sum(all_samples %in% existing_ids), "/", length(all_samples), "\n")
    cat("Missing objects: ", length(missing_ids), "\n")

    if (length(missing_ids) > 0) {
      cat("Samples to create:\n")
      print(missing_ids)
    } else {
      cat("All requested samples already have CellChat objects.\n")
    }
  }
}

# Run all requested samples.
# run_one() automatically skips an existing object and creates a
# missing object. Errors are logged without stopping other samples.
for (model_name in names(models)) {
  model_cfg <- models[[model_name]]

  for (group_definition in c("subtype4", "broadtype")) {

    cat("\n============================================================\n")
    cat("RUNNING:", model_name, "|", group_definition, "\n")
    cat("============================================================\n")

    for (sid in all_samples) {

      tryCatch(
        run_one(
          sample_id = sid,
          group_definition = group_definition,
          model_name = model_name,
          model_cfg = model_cfg
        ),
        error = function(e) {

          error_file <- file.path(
            log_out_dir,
            paste0(
              "ERROR_",
              model_name,
              "_",
              group_definition,
              "_",
              sid,
              ".txt"
            )
          )

          writeLines(
            c(
              paste0("Sample: ", sid),
              paste0("Model: ", model_name),
              paste0("Group: ", group_definition),
              paste0("Error: ", conditionMessage(e))
            ),
            error_file
          )

          message(
            "ERROR: ",
            model_name,
            " / ",
            group_definition,
            " / ",
            sid,
            " : ",
            conditionMessage(e)
          )
        }
      )
    }
  }
}

# ------------------------------------------------------------
# FINAL VERIFICATION
# ------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL CELLCHAT OBJECT VERIFICATION\n")
cat("============================================================\n")

for (model_name in names(models)) {
  for (group_definition in c("subtype4", "broadtype")) {

    obj_dir <- file.path(
      cellchat_out_dir,
      model_name,
      group_definition
    )

    existing_files <- list.files(
      obj_dir,
      pattern = "^CCI_CellChat_.*\\.rds$",
      full.names = FALSE
    )

    existing_ids <- sub(
      paste0("^CCI_CellChat_", group_definition, "_"),
      "",
      existing_files
    )
    existing_ids <- sub("\\.rds$", "", existing_ids)

    still_missing <- setdiff(all_samples, existing_ids)

    healthy_objects <- intersect(
      healthy_samples,
      existing_ids
    )

    risk_objects <- intersect(
      risk_samples,
      existing_ids
    )

    cat("\nModel:", model_name, "| Group:", group_definition, "\n")
    cat("Healthy objects:", length(healthy_objects), "/", length(healthy_samples), "\n")
    cat("Risk objects:   ", length(risk_objects), "/", length(risk_samples), "\n")
    cat("Total requested:", length(all_samples), "\n")
    cat("Still missing:  ", length(still_missing), "\n")

    if (length(still_missing) > 0) {
      cat("Still missing samples:\n")
      print(still_missing)
    } else {
      cat("PASS: All requested samples have CellChat objects.\n")
    }
  }
}

# ------------------------------------------------------------
# CONFIGURATION / SESSION
# ------------------------------------------------------------

config <- list(
  h5ad_file = h5ad_file,
  meta_file = meta_file,
  healthy_fp = healthy_fp,
  risk_fp = risk_fp,
  sample_column = SAMPLE_COL,
  broad_column = BROAD_COL,
  subtype_column = SUBTYPE_COL,
  expression_layer = EXPRESSION_LAYER,
  min_cells = MIN_CELLS,
  primary_model = list(type = PRIMARY_TYPE, trim = PRIMARY_TRIM),
  run_trimean_sensitivity = RUN_TRIMEAN_SENSITIVITY,
  population_size = POPULATION_SIZE,
  raw_use = RAW_USE,
  nboot = NBOOT,
  seed = SEED,
  samples_analyzed = all_samples
)

saveRDS(config, file.path(config_out_dir, "CellChat_run_config.rds"))
writeLines(capture.output(sessionInfo()), file.path(config_out_dir, "sessionInfo.txt"))

cat("\n============================================================\n")
cat("CELLCHAT COMPLETE\n")
cat("Output:", normalizePath(base_out_dir, mustWork = FALSE), "\n")
cat("Healthy files:", length(healthy_samples), "| Risk files:", length(risk_samples), "\n")
cat("Samples present in h5ad:", length(all_samples), "\n")
cat("============================================================\n")
