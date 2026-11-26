#!/usr/bin/env Rscript

# ============================================================
# 02_differential_CCI_AHI1_FINAL.R
#
# Small-cohort differential CellChat analysis.
#
# PRIMARY ANALYSIS
#   - truncatedMean10 CellChat objects
#   - subtype4 by default
#   - prenatal by default
#   - AHI1 risk vs healthy
#   - Model formula includes AHI1 status only.
#
# IMPORTANT FOR SMALL N
#   - No requirement for >10 samples per LR/pathway.
#   - Requires >=3 measured samples per AHI1 group per setting.
#   - Missing/not-measured combinations remain NA, not zero.
#   - Saves detection/presence summaries separately from inferential results.
#   - Raw P values, BH FDR, effect size and coverage are all retained.
# ============================================================

suppressPackageStartupMessages({
  library(reshape2)
  library(dplyr)
  library(limma)
  library(data.table)
})

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

BASE_OUT_DIR <- "cellchat_results_01"
MODEL_NAME <- "truncatedMean10"
GROUP_DEFINITION <- "subtype4"
RESTRICT_STAGE <- "prenatal"   # use NULL for all stages

healthy_fp <- "AHI1_rs13197384_healthy_samples.txt"
risk_fp    <- "AHI1_rs13197384_risk_homo_samples.txt"
subject_tsv <- "brainvar2_single_cell_detailed_metadata_v1.3.tsv"

STAGE_COL <- "Period"

MIN_NON_NA_PER_GROUP <- 3L
MIN_LR_SAMPLES <- 4L
MIN_PATHWAY_SAMPLES <- 4L
MIN_PROBABILITY <- 0.001

# ------------------------------------------------------------
# HELPERS
# ------------------------------------------------------------

safe_write_csv <- function(x, file) {
  data.table::fwrite(as.data.frame(x), file = file)
}

array_to_long <- function(arr, feature_name) {
  if (is.null(arr)) return(data.frame())
  if (length(dim(arr)) != 3L) {
    stop("Expected a 3D array, got: ", paste(dim(arr), collapse = " x "))
  }
  dn <- dimnames(arr)
  names(dn) <- c("Sender", "Receiver", feature_name)
  dimnames(arr) <- dn
  as.data.frame.table(arr, responseName = "value", stringsAsFactors = FALSE)
}

choose_covariates <- function(coldata, requested) {
  covars <- requested[requested %in% colnames(coldata)]
  if (length(covars) == 0) return(character(0))
  covars[sapply(coldata[covars], function(x) length(unique(na.omit(x))) > 1)]
}

fit_limma_smalln <- function(mat, coldata, prefix, out_dir) {
  if (is.null(mat) || nrow(mat) == 0) return(NULL)

  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  coverage <- t(vapply(seq_len(nrow(mat)), function(i) {
    present <- !is.na(mat[i, ])
    c(
      healthy = sum(present & coldata$AHI1_status == "healthy"),
      risk = sum(present & coldata$AHI1_status == "risk")
    )
  }, numeric(2)))

  coverage <- as.data.frame(coverage)
  coverage$settings <- rownames(mat)
  coverage$healthy_ok <- coverage$healthy >= MIN_NON_NA_PER_GROUP
  coverage$risk_ok <- coverage$risk >= MIN_NON_NA_PER_GROUP
  coverage$coverage_ok <- coverage$healthy_ok & coverage$risk_ok

  safe_write_csv(
    coverage,
    file.path(out_dir, paste0(prefix, "_coverage.csv"))
  )

  keep <- coverage$coverage_ok
  mat <- mat[keep, , drop = FALSE]
  if (nrow(mat) == 0) return(NULL)

  # Bounded communication probabilities -> arcsine-square-root transform.
  mat_trans <- asin(sqrt(pmax(pmin(mat, 1), 0)))

  safe_write_csv(
    data.frame(
      Sample = rownames(coldata),
      retained_for_model = TRUE,
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, paste0(prefix, "_complete_case_samples.csv"))
  )

  coldata$AHI1_status <- factor(
    coldata$AHI1_status,
    levels = c("healthy", "risk")
  )

  if (sum(coldata$AHI1_status == "healthy") < MIN_NON_NA_PER_GROUP ||
      sum(coldata$AHI1_status == "risk") < MIN_NON_NA_PER_GROUP) {
    warning(prefix, ": fewer than ", MIN_NON_NA_PER_GROUP, " complete samples per group.")
  }

  # Primary differential model: AHI1_status only.
  design <- model.matrix(~ AHI1_status, data = coldata)

  n_samples <- nrow(design)
  n_params <- ncol(design)
  residual_df <- n_samples - n_params

  safe_write_csv(
    data.frame(
      prefix = prefix,
      n_samples = n_samples,
      n_parameters = n_params,
      residual_df = residual_df,
      covariates = ifelse(length(covariates) == 0, "none", paste(covariates, collapse = ",")),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, paste0(prefix, "_model_info.csv"))
  )

  if (residual_df < 4) {
    warning(prefix, ": very low residual degrees of freedom: ", residual_df)
  }

  fit <- limma::lmFit(mat_trans, design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  coef_name <- grep("^AHI1_status", colnames(design), value = TRUE)[1]
  if (is.na(coef_name)) stop("Could not identify AHI1_status coefficient.")

  res <- limma::topTable(
    fit,
    coef = coef_name,
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  res$settings <- rownames(res)
  res$Model <- MODEL_NAME
  res$GroupDefinition <- GROUP_DEFINITION
  res
}

# ------------------------------------------------------------
# INPUT CHECKS
# ------------------------------------------------------------

for (f in c(healthy_fp, risk_fp, subject_tsv)) {
  if (!file.exists(f)) stop("Required file not found: ", f)
}

healthy_samples <- unique(trimws(readLines(healthy_fp)))
risk_samples <- unique(trimws(readLines(risk_fp)))
healthy_samples <- healthy_samples[nzchar(healthy_samples)]
risk_samples <- risk_samples[nzchar(risk_samples)]

if (length(intersect(healthy_samples, risk_samples)) > 0) {
  stop("Overlap between healthy and risk sample lists.")
}

cat("Healthy samples in list:", length(healthy_samples), "\n")
cat("Risk samples in list:", length(risk_samples), "\n")

# ------------------------------------------------------------
# SAMPLE METADATA
# Only sample ID and developmental stage are needed.
# ------------------------------------------------------------

subject_raw <- data.table::fread(subject_tsv)
subject_raw <- as.data.frame(subject_raw, stringsAsFactors = FALSE)

subject_meta <- subject_raw %>%
  transmute(
    Sample = as.character(.data[["BioSampleID"]]),
    stage = .data[[STAGE_COL]]
  ) %>%
  distinct(Sample, .keep_all = TRUE) %>%
  mutate(
    stage_num = suppressWarnings(as.numeric(as.character(stage))),
    stage_group = case_when(
      stage_num %in% 1:6 ~ "prenatal",
      stage_num == 7 ~ NA_character_,
      stage_num %in% 8:15 ~ "postnatal",
      TRUE ~ NA_character_
    )
  )

sample_meta <- data.frame(
  Sample = unique(c(healthy_samples, risk_samples)),
  stringsAsFactors = FALSE
) %>%
  mutate(
    AHI1_status = case_when(
      Sample %in% healthy_samples ~ "healthy",
      Sample %in% risk_samples ~ "risk",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(subject_meta, by = "Sample") %>%
  filter(!is.na(AHI1_status))

if (!is.null(RESTRICT_STAGE)) {
  sample_meta <- sample_meta %>%
    filter(stage_group == RESTRICT_STAGE)
}

# ------------------------------------------------------------
# LOAD CELLCHAT OBJECTS
# ------------------------------------------------------------

analysis_dir <- file.path(
  BASE_OUT_DIR,
  "cellchat_objects",
  MODEL_NAME,
  GROUP_DEFINITION
)

if (!dir.exists(analysis_dir)) {
  stop("CellChat directory not found: ", analysis_dir)
}

all_files <- list.files(
  analysis_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No CellChat RDS files found in: ", analysis_dir)
}

CCI_sample <- lapply(all_files, readRDS)

# Get sample IDs from object metadata rather than file-name parsing.
obj_sample_ids <- vapply(CCI_sample, function(x) {
  if (!is.null(x@meta$samples)) {
    s <- unique(as.character(x@meta$samples))
    if (length(s) == 1) return(s)
  }
  NA_character_
}, character(1))

keep_obj <- !is.na(obj_sample_ids)
CCI_sample <- CCI_sample[keep_obj]
obj_sample_ids <- obj_sample_ids[keep_obj]
names(CCI_sample) <- obj_sample_ids

keep_samples <- intersect(names(CCI_sample), sample_meta$Sample)
CCI_sample <- CCI_sample[keep_samples]
sample_meta <- sample_meta[match(keep_samples, sample_meta$Sample), , drop = FALSE]
rownames(sample_meta) <- sample_meta$Sample

if (length(CCI_sample) == 0) {
  stop("No CellChat objects matched the sample metadata after stage restriction.")
}

cat("\nFinal samples entering differential analysis:\n")
print(table(sample_meta$AHI1_status))

out_dir <- file.path(
  BASE_OUT_DIR,
  "homo_differential",
  MODEL_NAME,
  GROUP_DEFINITION,
  if (is.null(RESTRICT_STAGE)) "all_stages" else RESTRICT_STAGE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

safe_write_csv(sample_meta, file.path(out_dir, "sample_metadata_used.csv"))

# ------------------------------------------------------------
# EXTRACT LR / PATHWAY ARRAYS
# ------------------------------------------------------------

lr_list <- lapply(CCI_sample, function(x) x@net$prob)
pathway_list <- lapply(CCI_sample, function(x) {
  if (!is.null(x@netP$prob)) x@netP$prob else NULL
})

# Detected pathways by sample.
detected_pathways <- do.call(rbind, lapply(names(pathway_list), function(sid) {
  arr <- pathway_list[[sid]]
  if (is.null(arr)) return(NULL)
  data.frame(
    Sample = sid,
    Pathway = dimnames(arr)[[3]],
    stringsAsFactors = FALSE
  )
}))

if (is.null(detected_pathways)) {
  detected_pathways <- data.frame(Sample = character(), Pathway = character())
}

safe_write_csv(
  detected_pathways,
  file.path(out_dir, "pathways_detected_by_sample.csv")
)

# ------------------------------------------------------------
# LR LONG TABLE + DETECTION
# ------------------------------------------------------------

melt_lr <- function(arr, sid) {
  if (is.null(arr)) return(data.frame())
  df <- array_to_long(arr, "LR")
  df$Sample <- sid
  df
}

df_lr <- do.call(rbind, Map(melt_lr, lr_list, names(lr_list)))
rownames(df_lr) <- NULL

if (nrow(df_lr) == 0) {
  stop("No LR probability rows were extracted.")
}

df_lr$settings <- paste(df_lr$Sender, df_lr$Receiver, df_lr$LR, sep = "|")
safe_write_csv(df_lr, file.path(out_dir, "LR_all_samples_raw.csv"))

lr_detection <- df_lr %>%
  group_by(LR) %>%
  summarise(
    n_samples_nonzero = sum(!is.na(value) & value > 0),
    n_samples_above_threshold = sum(!is.na(value) & value > MIN_PROBABILITY),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples_above_threshold), desc(n_samples_nonzero))

safe_write_csv(lr_detection, file.path(out_dir, "LR_detection_summary_before_filter.csv"))

lr_diag <- lr_detection[
  grepl("Glutamate|Glu-|GRIA|GRIK|GRM|SLC17|SLC1A|GLS|GLS2|Acetylcholine|Ach-|CHAT|SLC5A7|SLC18A3|CHRN|CHRNA|CHRM",
        lr_detection$LR, ignore.case = TRUE),
  , drop = FALSE
]
safe_write_csv(lr_diag, file.path(out_dir, "LR_detection_glutamate_cholinergic.csv"))

keep_lr <- lr_detection$LR[
  lr_detection$n_samples_nonzero >= MIN_LR_SAMPLES
]

df_lr_filt <- df_lr[df_lr$LR %in% keep_lr, , drop = FALSE]

setting_detection <- df_lr_filt %>%
  group_by(settings) %>%
  summarise(
    n_samples_above_threshold = sum(!is.na(value) & value > MIN_PROBABILITY),
    .groups = "drop"
  ) %>%
  arrange(desc(n_samples_above_threshold))

safe_write_csv(
  setting_detection,
  file.path(out_dir, "LR_setting_detection_summary_before_filter.csv")
)

keep_settings <- setting_detection$settings[
  setting_detection$n_samples_above_threshold >= MIN_LR_SAMPLES
]

df_lr_filt <- df_lr_filt[
  df_lr_filt$settings %in% keep_settings,
  , drop = FALSE
]

# ------------------------------------------------------------
# LR MATRIX / DIFFERENTIAL TEST
# ------------------------------------------------------------

if (nrow(df_lr_filt) > 0) {
  cci_mat_df <- reshape2::dcast(
    df_lr_filt,
    settings ~ Sample,
    value.var = "value",
    fill = NA
  )
  rownames(cci_mat_df) <- cci_mat_df$settings
  cci_mat_df$settings <- NULL
  cci_mat <- as.matrix(cci_mat_df)
  cci_mat <- cci_mat[, sample_meta$Sample, drop = FALSE]
} else {
  cci_mat <- matrix(
    numeric(0),
    nrow = 0,
    ncol = nrow(sample_meta),
    dimnames = list(NULL, sample_meta$Sample)
  )
}

lr_results <- fit_limma_smalln(
  cci_mat,
  sample_meta,
  "LR",
  out_dir
)

if (!is.null(lr_results)) {
  safe_write_csv(
    lr_results,
    file.path(out_dir, "AHI1_differential_LR_settings.csv")
  )

  target <- lr_results[
    grepl("Glutamate|Glu-|GRIA|GRIK|GRM|SLC17|SLC1A|GLS|GLS2|Acetylcholine|Ach-|CHAT|SLC5A7|SLC18A3|CHRN|CHRNA|CHRM",
          lr_results$settings, ignore.case = TRUE),
    , drop = FALSE
  ]
  safe_write_csv(
    target,
    file.path(out_dir, "AHI1_differential_glutamate_cholinergic_LR.csv")
  )
}

# ------------------------------------------------------------
# PATHWAY LONG TABLE + DETECTION
# ------------------------------------------------------------

melt_pathway <- function(arr, sid) {
  if (is.null(arr)) return(data.frame())
  df <- array_to_long(arr, "Pathway")
  df$Sample <- sid
  df
}

df_path <- do.call(rbind, Map(melt_pathway, pathway_list, names(pathway_list)))
rownames(df_path) <- NULL

path_results <- NULL
path_global_results <- NULL
path_mat <- matrix(
  numeric(0),
  nrow = 0,
  ncol = nrow(sample_meta),
  dimnames = list(NULL, sample_meta$Sample)
)

if (nrow(df_path) > 0) {
  df_path$settings <- paste(df_path$Sender, df_path$Receiver, df_path$Pathway, sep = "|")
  safe_write_csv(df_path, file.path(out_dir, "pathway_all_samples_raw.csv"))

  pathway_detection <- df_path %>%
    group_by(Pathway) %>%
    summarise(
      n_samples_nonzero = sum(!is.na(value) & value > 0),
      n_samples_above_threshold = sum(!is.na(value) & value > MIN_PROBABILITY),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples_above_threshold), desc(n_samples_nonzero))

  safe_write_csv(
    pathway_detection,
    file.path(out_dir, "pathway_detection_summary_before_filter.csv")
  )

  pathway_diag <- pathway_detection[
    grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline",
          pathway_detection$Pathway, ignore.case = TRUE),
    , drop = FALSE
  ]
  safe_write_csv(
    pathway_diag,
    file.path(out_dir, "pathway_detection_glutamate_cholinergic.csv")
  )

  path_setting_detection <- df_path %>%
    group_by(settings) %>%
    summarise(
      n_samples_above_threshold = sum(!is.na(value) & value > MIN_PROBABILITY),
      .groups = "drop"
    ) %>%
    arrange(desc(n_samples_above_threshold))

  safe_write_csv(
    path_setting_detection,
    file.path(out_dir, "pathway_setting_detection_summary_before_filter.csv")
  )

  keep_path_settings <- path_setting_detection$settings[
    path_setting_detection$n_samples_above_threshold >= MIN_PATHWAY_SAMPLES
  ]

  df_path_filt <- df_path[
    df_path$settings %in% keep_path_settings,
    , drop = FALSE
  ]

  if (nrow(df_path_filt) > 0) {
    path_mat_df <- reshape2::dcast(
      df_path_filt,
      settings ~ Sample,
      value.var = "value",
      fill = NA
    )
    rownames(path_mat_df) <- path_mat_df$settings
    path_mat_df$settings <- NULL
    path_mat <- as.matrix(path_mat_df)
    path_mat <- path_mat[, sample_meta$Sample, drop = FALSE]

    path_results <- fit_limma_smalln(
      path_mat,
      sample_meta,
      "pathway_sender_receiver",
      out_dir
    )

    if (!is.null(path_results)) {
      safe_write_csv(
        path_results,
        file.path(out_dir, "AHI1_differential_pathway_sender_receiver.csv")
      )

      target_path <- path_results[
        grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline",
              path_results$settings, ignore.case = TRUE),
        , drop = FALSE
      ]
      safe_write_csv(
        target_path,
        file.path(out_dir, "AHI1_differential_glutamate_cholinergic_pathway_sender_receiver.csv")
      )
    }
  }

  # Global pathway totals per sample.
  global_df <- df_path %>%
    group_by(Pathway, Sample) %>%
    summarise(total_probability = sum(value, na.rm = TRUE), .groups = "drop") %>%
    reshape2::dcast(Pathway ~ Sample, value.var = "total_probability", fill = 0)

  global_mat <- as.matrix(global_df[, -1, drop = FALSE])
  rownames(global_mat) <- global_df$Pathway
  global_mat <- global_mat[, sample_meta$Sample, drop = FALSE]

  path_global_results <- fit_limma_smalln(
    global_mat,
    sample_meta,
    "pathway_global",
    out_dir
  )

  if (!is.null(path_global_results)) {
    safe_write_csv(
      path_global_results,
      file.path(out_dir, "AHI1_differential_pathway_global.csv")
    )

    target_global <- path_global_results[
      grepl("^Glutamate$|^Ach$|acetylcholine|cholinergic|choline",
            path_global_results$settings, ignore.case = TRUE),
      , drop = FALSE
    ]
    safe_write_csv(
      target_global,
      file.path(out_dir, "AHI1_differential_glutamate_cholinergic_pathway_global.csv")
    )
  }
}

# ------------------------------------------------------------
# SAVE PLOT INPUTS / CONFIG
# ------------------------------------------------------------

saveRDS(
  list(
    cci_mat = cci_mat,
    path_mat = path_mat,
    coldata = sample_meta,
    lr_results = lr_results,
    pathway_results = path_results,
    pathway_global_results = path_global_results,
    model = MODEL_NAME,
    group_definition = GROUP_DEFINITION,
    restrict_stage = RESTRICT_STAGE,
    min_non_na_per_group = MIN_NON_NA_PER_GROUP,
    min_lr_samples = MIN_LR_SAMPLES,
    min_pathway_samples = MIN_PATHWAY_SAMPLES,
    min_probability = MIN_PROBABILITY
  ),
  file.path(out_dir, "AHI1_plot_inputs.rds")
)

safe_write_csv(
  data.frame(
    parameter = c(
      "MODEL_NAME",
      "GROUP_DEFINITION",
      "RESTRICT_STAGE",
      "MIN_NON_NA_PER_GROUP",
      "MIN_LR_SAMPLES",
      "MIN_PATHWAY_SAMPLES",
      "MIN_PROBABILITY",
      "MODEL_TERMS"
    ),
    value = c(
      MODEL_NAME,
      GROUP_DEFINITION,
      ifelse(is.null(RESTRICT_STAGE), "NULL", RESTRICT_STAGE),
      MIN_NON_NA_PER_GROUP,
      MIN_LR_SAMPLES,
      MIN_PATHWAY_SAMPLES,
      MIN_PROBABILITY,
      "AHI1_status"
    ),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "analysis_configuration.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("DIFFERENTIAL CELLCHAT ANALYSIS COMPLETE\n")
cat("Output:", normalizePath(out_dir, mustWork = FALSE), "\n")
cat("============================================================\n")
