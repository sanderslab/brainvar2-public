#!/usr/bin/env Rscript

# =============================================================================
# GRN MODEL HELPERS
# =============================================================================
# This script performs GRN inference using X, Y, Z matrices and motif matching
# results from period-specific data retrieval. It follows Step 5 from
# complete_grn_pipeline.R but works with period-level aggregated data.
#
# Internal module shared by the phase and locus entry points.

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(Matrix)
  library(grpregOverlap)
  library(progress)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Paths are supplied by the public entry point when this module is sourced.
if (!exists("BASE_DIR", inherits = FALSE)) BASE_DIR <- getwd()
if (!exists("DEFAULT_INPUT_DIR", inherits = FALSE)) DEFAULT_INPUT_DIR <- file.path(BASE_DIR, "period_data_retrieval_results")
if (!exists("DEFAULT_OUTPUT_DIR", inherits = FALSE)) DEFAULT_OUTPUT_DIR <- file.path(BASE_DIR, "period_grn_inference_results")
if (!exists("PERIOD_MOTIF_RESULTS_DIR", inherits = FALSE)) PERIOD_MOTIF_RESULTS_DIR <- file.path(BASE_DIR, "period_motif_matching_results")

# =============================================================================
# GRN INFERENCE FUNCTIONS (same as batch version)
# =============================================================================

#' Construct groups based on motif binding matrix
#' @param B Binary motif binding matrix (TFs x CREs)
#' @return List of groups for group lasso
#' @note Only TF-CRE pairs are included in groups. Single TFs or CREs without bindings are ignored.
construct_group <- function(B){
  N_TF <- nrow(B)
  N_CRE <- ncol(B)

  group <- list()

  # Only create groups for TF-CRE pairs (where motif binding exists)
  # Single TFs or CREs without bindings are ignored
  binding_indices <- which(B == 1, arr.ind = TRUE)

  for (i in seq_len(nrow(binding_indices))){
    grp_name <- paste0("grp", i)
    tf_index <- binding_indices[i, 1]
    cre_index <- binding_indices[i, 2] + N_TF

    if (tf_index <= N_TF && cre_index <= N_TF + N_CRE) {
      group[[grp_name]] <- c(tf_index, cre_index)
    }
  }

  return(group)
}

#' Fit group lasso model
#' @param TF_RNA TF expression matrix (samples × TFs)
#' @param CRE_ATAC CRE accessibility matrix (samples × CREs)
#' @param TG_RNA Target gene expression vector (samples)
#' @param group Group structure
#' @param lambda Lambda parameter
#' @param verbose Verbosity level
#' @param seed Random seed
#' @return List with inference model and coefficients
grpreg_fit <- function(TF_RNA, CRE_ATAC, TG_RNA, group, batch_covariates = NULL,
                       lambda = NULL, verbose = 0, seed = 1234){
  # Ensure inputs are matrices
  TF_RNA <- as.matrix(TF_RNA)
  CRE_ATAC <- as.matrix(CRE_ATAC)

  # Combine TF and CRE matrices directly (both are samples × features)
  X_tf_cre <- cbind(TF_RNA, CRE_ATAC)
  X_tf_cre <- as.matrix(X_tf_cre)

  # Store number of TF/CRE variables (before adding batch covariates)
  n_tf_cre_vars <- ncol(X_tf_cre)
  n_tf_vars <- ncol(TF_RNA)
  n_cre_vars <- ncol(CRE_ATAC)

  # Add batch covariates if provided (these will be unpenalized - not in any group)
  if (!is.null(batch_covariates)) {
    batch_covariates <- as.matrix(batch_covariates)
    # Ensure batch covariates have same number of rows as X
    if (nrow(batch_covariates) != nrow(X_tf_cre)) {
      stop("ERROR: Batch covariates must have same number of rows as TF/CRE matrices")
    }
    X <- cbind(X_tf_cre, batch_covariates)
    if (verbose > 0) {
      cat("  Added", ncol(batch_covariates), "batch covariates (unpenalized)\n")
    }
  } else {
    X <- X_tf_cre
  }

  X <- as.matrix(X)

  # Check dimensions and validate group indices
  n_vars <- ncol(X)

  # Validate group indices don't exceed TF/CRE matrix dimensions (not including batch covariates)
  # Batch covariates are not in any group (unpenalized)
  if (length(group) > 0) {
    max_group_index <- max(unlist(group))
    if (max_group_index > n_tf_cre_vars) {
      warning("Group indices exceed TF/CRE matrix dimensions. Adjusting groups.")
      group <- lapply(group, function(g) {
        g[g <= n_tf_cre_vars]
      })
      group <- group[sapply(group, length) > 0]
    }

    # Ensure all group indices are within TF/CRE range (batch covariates come after)
    all_group_indices <- unique(unlist(group))
    if (any(all_group_indices > n_tf_cre_vars)) {
      stop("ERROR: Group indices must be <= number of TF/CRE variables (", n_tf_cre_vars, ")")
    }
  }

  # Normalize columns of X, handling zero variance
  # NOTE: Input data (X, Y, Z) should already be normalized/transformed for single-cell data:
  #   - RNA expression (X, Y): Typically log-normalized (log1p) or log2-transformed
  #   - ATAC accessibility (Z): May be normalized or raw counts depending on data source
  #   - If batch correction was applied, data may have been normalized via Seurat's NormalizeData (log normalization)
  # This standardization step (z-score) is for regression modeling, not data normalization
  X_scaled <- X
  zero_var_cols <- integer(0)
  for (j in seq_len(ncol(X))) {
    var_j <- var(X[,j], na.rm = TRUE)
    if (is.na(var_j) || var_j == 0 || !is.finite(var_j)) {
      zero_var_cols <- c(zero_var_cols, j)
      X_scaled[,j] <- X[,j]
    } else {
      X_scaled[,j] <- scale(X[,j])[,1]
    }
  }
  X <- X_scaled

  # Remove constant columns from group structure if needed
  if (length(zero_var_cols) > 0) {
    group_filtered <- lapply(group, function(g) {
      g[!g %in% zero_var_cols]
    })
    group_filtered <- group_filtered[sapply(group_filtered, length) > 0]
    group <- group_filtered
  }

  Y <- as.numeric(TG_RNA)

  # Check for valid data
  if (any(!is.finite(Y)) || any(!is.finite(X))) {
    stop("ERROR: Non-finite values in X or Y")
  }

  # Sample splitting for model selection and inference
  set.seed(seed)
  fold_indices <- sample(1: 2, size = nrow(X), replace = TRUE)

  X_sel <- X[fold_indices == 1, , drop = FALSE]
  Y_sel <- Y[fold_indices == 1]

  X_inf <- X[fold_indices == 2, , drop = FALSE]
  Y_inf <- Y[fold_indices == 2]

  # Ensure we have enough samples in each fold
  if (nrow(X_sel) < 5 || nrow(X_inf) < 5) {
    warning("Too few samples after splitting. Using all data for selection.")
    X_sel <- X
    Y_sel <- Y
    X_inf <- X
    Y_inf <- Y
  }

  # Validate group structure
  if (length(group) == 0) {
    stop("ERROR: Empty group structure after filtering")
  }

  # Check that all group indices are valid (only for TF/CRE variables, not batch covariates)
  all_group_indices <- unique(unlist(group))
  if (any(all_group_indices < 1) || any(all_group_indices > n_tf_cre_vars)) {
    stop("ERROR: Invalid group indices. Must be between 1 and", n_tf_cre_vars)
  }

  # Tune the model with error handling
  lambda_seq <- NULL
  tryCatch({
    tune <- cv.grpregOverlap(X_sel, Y_sel, group, penalty = 'grLasso')
  }, error = function(e) {
    cat("  Warning: cv.grpregOverlap failed, computing lambda sequence manually...\n")
    tryCatch({
      fit_temp <- grpregOverlap(X_sel, Y_sel, group, penalty = 'grLasso',
                               nlambda = 50, lambda.min = 0.001)
      lambda_seq <<- fit_temp$lambda
      cat("  Computed lambda sequence:", length(lambda_seq), "values\n")
    }, error = function(e2) {
      cat("  Warning: Automatic lambda computation failed, using manual sequence...\n")
      lambda_seq <<- exp(seq(log(0.1), log(0.001), length.out = 20))
      cat("  Using manual lambda sequence:", length(lambda_seq), "values\n")
    })

    if (!is.null(lambda_seq)) {
      tryCatch({
        tune <<- cv.grpregOverlap(X_sel, Y_sel, group, penalty = 'grLasso',
                                 lambda = lambda_seq)
      }, error = function(e3) {
        cat("  Warning: cv.grpregOverlap still failed, using grpregOverlap without CV...\n")
        if (is.null(lambda_seq)) {
          lambda_seq <<- exp(seq(log(0.1), log(0.001), length.out = 20))
        }
        tune <<- grpregOverlap(X_sel, Y_sel, group, penalty = 'grLasso',
                              lambda = lambda_seq)
      })
    } else {
      cat("  Warning: Using grpregOverlap without CV as last resort...\n")
      lambda_seq <<- exp(seq(log(0.1), log(0.001), length.out = 20))
      tune <<- grpregOverlap(X_sel, Y_sel, group, penalty = 'grLasso',
                            lambda = lambda_seq)
    }
  })

  coef_sel <- tryCatch({
    coef(tune)
  }, error = function(e) {
    if (inherits(tune, "cv.grpregOverlap")) {
      coef(tune, lambda = tune$lambda.min)
    } else {
      coef(tune)
    }
  })

  coef_sel_latent <- tryCatch({
    coef(tune, latent = TRUE)
  }, error = function(e) {
    coef_sel
  })

  # Extract selected variables (excluding intercept)
  if (length(coef_sel) > 1) {
    var_sel <- which(coef_sel[2:length(coef_sel)] != 0)
  } else {
    var_sel <- integer(0)
  }

  # Perform inference on the other fold
  # Note: var_sel indices refer only to TF/CRE variables (first n_tf_cre_vars columns)
  # Batch covariates are in columns (n_tf_cre_vars + 1):ncol(X) and should always be included

  if (!is.null(batch_covariates)) {
    # Get batch covariates for inference fold
    batch_cov_inf <- batch_covariates[fold_indices == 2, , drop = FALSE]
    if (nrow(batch_cov_inf) < 5) {
      # If too few samples in fold 2, use all data
      batch_cov_inf <- batch_covariates
    }

    if (length(var_sel) > 0) {
      # Select TF/CRE variables (var_sel indices are 1-based for TF/CRE only)
      X_inf_tf_cre <- X_inf[, var_sel, drop = FALSE]
      # Always include batch covariates (unpenalized)
      X_inf_sel <- cbind(X_inf_tf_cre, batch_cov_inf)
      df <- data.frame(Y_inf, X_inf_sel)
      inf_model <- lm(Y_inf ~ ., data = df)
    } else {
      # Only batch covariates (no TF/CRE selected)
      df <- data.frame(Y_inf, batch_cov_inf)
      inf_model <- lm(Y_inf ~ ., data = df)
    }
  } else {
    # No batch covariates
    if (length(var_sel) > 0) {
      X_inf_sel <- X_inf[, var_sel, drop = FALSE]
      df <- data.frame(Y_inf, X_inf_sel)
      inf_model <- lm(Y_inf ~ ., data = df)
    } else {
      inf_model <- lm(Y_inf ~ 1)
    }
  }

  return(list(
    inference_model = inf_model,
    selection_groups_coef = coef_sel_latent,
    n_selected_variables = length(var_sel)
  ))
}

#' Find significant TFs
#' @param combined_p_vals Combined p-values (should be BH-adjusted)
#' @param motif_matrix Motif binding matrix
#' @param p_value_cutoff P-value cutoff (for adjusted p-values)
#' @return Vector of significant TFs
find_significant_tfs <- function(combined_p_vals, motif_matrix, p_value_cutoff = 0.05){
  cat("  Identifying significant TFs from adjusted p-values...\n")
  cat("  Total variables with p-values:", length(combined_p_vals), "\n")
  cat("  Variables with adjusted p <", p_value_cutoff, ":", sum(combined_p_vals < p_value_cutoff, na.rm = TRUE), "\n")

  sig_combined_p_vals <- combined_p_vals[combined_p_vals < p_value_cutoff & !is.na(combined_p_vals)]

  if (length(sig_combined_p_vals) == 0) {
    cat("  WARNING: No variables have significant p-values!\n")
    cat("  Minimum p-value:", min(combined_p_vals), "\n")
    cat("  Median p-value:", median(combined_p_vals), "\n")
    return(character(0))
  }

  # Extract TF and CRE names only from the modeled biological features.
  # Unpenalized covariates such as period_* or batch_* must not be treated as CREs.
  all_variable_names <- names(combined_p_vals)
  motif_tfs <- rownames(motif_matrix)
  motif_cres <- colnames(motif_matrix)
  potential_tfs <- intersect(all_variable_names, motif_tfs)
  potential_cres <- intersect(all_variable_names, motif_cres)

  cat("  Potential TFs in variable names:", length(potential_tfs), "\n")
  cat("  Potential CREs in variable names:", length(potential_cres), "\n")

  # Strategy: TFs must be significant AND bind to a significant CRE
  # A TF is identified as significant only if:
  #   1. The TF itself has a significant p-value (adjusted_p_value < cutoff)
  #   2. The TF binds to a CRE that also has a significant p-value
  sig_cre_names <- potential_cres[potential_cres %in% names(sig_combined_p_vals)]
  sig_cre_names_clean <- gsub("\\.", "-", sig_cre_names)

  sig_tf_via_cre <- character(0)
  if (length(sig_cre_names_clean) > 0) {
    matched_cre_cols <- match(sig_cre_names_clean, colnames(motif_matrix))
    matched_cre_cols <- matched_cre_cols[!is.na(matched_cre_cols)]

    if (length(matched_cre_cols) > 0) {
      for (cre_idx in matched_cre_cols) {
        cre_name <- colnames(motif_matrix)[cre_idx]
        binding_tfs <- rownames(motif_matrix)[motif_matrix[, cre_idx] == 1]
        # Only include TFs that are also significant (have significant p-values)
        binding_tfs_sig <- binding_tfs[binding_tfs %in% names(sig_combined_p_vals)]
        sig_tf_via_cre <- c(sig_tf_via_cre, binding_tfs_sig)
      }
    } else {
      # Try fuzzy matching for CRE names
      for (cre_name in sig_cre_names_clean) {
        cre_matches <- which(colnames(motif_matrix) == cre_name)
        if (length(cre_matches) == 0) {
          cre_parts <- strsplit(cre_name, "-")[[1]]
          if (length(cre_parts) >= 3) {
            cre_matches <- which(
              startsWith(colnames(motif_matrix), paste0(cre_parts[1], "-")) &
              grepl(cre_parts[2], colnames(motif_matrix)) &
              grepl(cre_parts[3], colnames(motif_matrix))
            )
          }
        }
        if (length(cre_matches) > 0) {
          for (cre_idx in cre_matches) {
            binding_tfs <- rownames(motif_matrix)[motif_matrix[, cre_idx] == 1]
            # Only include TFs that are also significant (have significant p-values)
            binding_tfs_sig <- binding_tfs[binding_tfs %in% names(sig_combined_p_vals)]
            sig_tf_via_cre <- c(sig_tf_via_cre, binding_tfs_sig)
          }
        }
      }
    }
  }

  # Final significant TFs: only those that bind to significant CREs AND are themselves significant
  final_sig_tfs <- unique(sig_tf_via_cre)
  final_sig_tfs <- sort(final_sig_tfs)

  cat("  Significant TFs (via CRE binding, TF must also be significant):", length(final_sig_tfs), "\n")
  if (length(final_sig_tfs) > 0) {
    cat("    ", paste(final_sig_tfs, collapse = ", "), "\n")
  }

  return(final_sig_tfs)
}

#' Cauchy combination test
#' @param p_values Vector of p-values
#' @return Combined p-value
cauchy_combine_p <- function(p_values){
  # Remove NA values before combining
  p_values_clean <- p_values[!is.na(p_values)]

  # If all values are NA, return NA
  if (length(p_values_clean) == 0) {
    return(NA_real_)
  }

  # If some values are NA, use only the non-NA values for combination
  # This is reasonable because NA typically indicates the variable was not
  # selected or the inference model couldn't compute a p-value (e.g., perfect collinearity)
  cauchy_transformed_p <- tan((0.5 - p_values_clean) * pi)
  test_statistics <- sum(cauchy_transformed_p) / length(p_values_clean)
  combined_p_value <- 1/2 - atan(test_statistics) / pi
  return(combined_p_value)
}

#' Run GRN inference
#' @param X_matrix TF expression matrix
#' @param Y_matrix Target gene expression matrix
#' @param Z_matrix CRE accessibility matrix
#' @param motif_matrix Motif binding matrix
#' @param batch_covariates Optional batch dummy variables matrix (samples × batches) - unpenalized
#' @param n_runs Number of bootstrap runs
#' @param p_value_cutoff P-value cutoff
#' @return List with significant TFs and results
classify_grn_variable <- function(variables, motif_matrix) {
  ifelse(
    variables %in% rownames(motif_matrix),
    "TF",
    ifelse(grepl("^period_", variables), "covariate", "CRE")
  )
}

safe_selected_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_selected_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

run_grn_inference <- function(X_matrix, Y_matrix, Z_matrix, motif_matrix,
                              batch_covariates = NULL,
                              n_runs = 100, p_value_cutoff = 0.05,
                              progress_file = NULL,
                              progress_every = 5) {
  cat("\nRunning GRN inference...\n")

  # Prepare data (matrices are already in samples × features format)
  TG_RNA <- as.numeric(Y_matrix[,1])

  # Construct groups
  group <- construct_group(motif_matrix)
  cat("✓ Constructed", length(group), "groups for group lasso\n")

  if (!is.null(batch_covariates)) {
    cat("  Batch covariates will be included as unpenalized variables\n")
  }

  # Run bootstrap iterations
  cat("Running", n_runs, "bootstrap iterations...\n")
  pb <- progress_bar$new(total = n_runs)
  if (!is.null(progress_file)) {
    writeLines(sprintf("0/%d", n_runs), progress_file)
  }

  p_val_list <- list()
  run_seconds <- numeric(n_runs)
  run_selected_vars <- integer(n_runs)
  for (seed in seq_len(n_runs)) {
    run_t0 <- proc.time()[["elapsed"]]
    fit <- grpreg_fit(X_matrix, Z_matrix, TG_RNA, group, batch_covariates = batch_covariates,
                     verbose = 0, seed = seed)
    run_t1 <- proc.time()[["elapsed"]]
    run_seconds[seed] <- run_t1 - run_t0
    run_selected_vars[seed] <- ifelse(is.null(fit$n_selected_variables), NA_integer_, fit$n_selected_variables)

    inference_model <- fit$inference_model
    coef_summary <- summary(inference_model)$coefficients

    # Extract post-selection statistics, excluding the intercept.
    if (nrow(coef_summary) > 1) {
      p_values <- coef_summary[2:nrow(coef_summary), 4]  # Remove intercept
      coefficients <- coef_summary[2:nrow(coef_summary), 1]
    } else {
      # Only intercept in model (no variables selected)
      p_values <- numeric(0)
      coefficients <- numeric(0)
      names(p_values) <- character(0)
      names(coefficients) <- character(0)
    }

    # FIX: Normalize variable names (R's data.frame converts dashes to dots)
    # Convert dots back to dashes for CRE names to match motif matrix
    var_names <- names(p_values)
    if (length(var_names) > 0) {
      var_names_normalized <- gsub("\\.", "-", var_names)
      names(p_values) <- var_names_normalized
      names(coefficients) <- var_names_normalized
    } else {
      var_names_normalized <- character(0)
    }

    p_val_list[[seed]] <- list(
      variables = var_names_normalized,
      p_values = p_values,
      coefficients = coefficients,
      elapsed_seconds = run_seconds[seed],
      n_selected_variables = run_selected_vars[seed]
    )

    pb$tick()
    if ((seed %% progress_every) == 0 || seed == 1 || seed == n_runs) {
      cat(sprintf("  Bootstrap progress: %d/%d\n", seed, n_runs))
    }
    cat(sprintf("    run=%d elapsed=%.3fs selected_vars=%d\n", seed, run_seconds[seed], run_selected_vars[seed]))
    if (!is.null(progress_file)) {
      writeLines(sprintf("%d/%d", seed, n_runs), progress_file)
    }
  }
  cat(sprintf("  Bootstrap timing summary (seconds): min=%.3f median=%.3f mean=%.3f max=%.3f\n",
              min(run_seconds), median(run_seconds), mean(run_seconds), max(run_seconds)))
  cat(sprintf("  Selected-variable summary: min=%d median=%.1f mean=%.2f max=%d\n",
              min(run_selected_vars), median(run_selected_vars), mean(run_selected_vars), max(run_selected_vars)))

  # Combine p-values
  all_variables <- unique(unlist(lapply(p_val_list, function(x) x$variables)))
  if (length(all_variables) == 0) {
    stop("No variables were selected in any bootstrap run")
  }
  p_val_matrix <- matrix(0.95, nrow = length(all_variables), ncol = n_runs)
  coefficient_matrix <- matrix(NA_real_, nrow = length(all_variables), ncol = n_runs)
  rownames(p_val_matrix) <- all_variables
  rownames(coefficient_matrix) <- all_variables

  for (seed in seq_len(n_runs)) {
    run_data <- p_val_list[[seed]]
    for (var in run_data$variables) {
      p_val_matrix[var, seed] <- run_data$p_values[var]
      coefficient_matrix[var, seed] <- run_data$coefficients[var]
    }
  }

  combined_p_values <- apply(p_val_matrix, 1, cauchy_combine_p)

  # Apply Benjamini-Hochberg (BH) FDR correction for multiple testing
  n_tests <- length(combined_p_values)
  adjusted_p_values <- p.adjust(combined_p_values, method = "BH")
  cat("  Applied Benjamini-Hochberg (BH) FDR correction for", n_tests, "tests\n")

  bootstrap_details <- do.call(rbind, lapply(seq_len(n_runs), function(seed) {
    run_data <- p_val_list[[seed]]
    if (length(run_data$variables) == 0) {
      return(data.frame(
        run = integer(0), variable = character(0), variable_type = character(0),
        coefficient = numeric(0), post_selection_p_value = numeric(0),
        stringsAsFactors = FALSE
      ))
    }
    data.frame(
      run = seed,
      variable = run_data$variables,
      variable_type = classify_grn_variable(run_data$variables, motif_matrix),
      coefficient = as.numeric(run_data$coefficients),
      post_selection_p_value = as.numeric(run_data$p_values),
      stringsAsFactors = FALSE
    )
  }))

  score_summary <- data.frame(
    Variable = all_variables,
    Variable_Type = classify_grn_variable(all_variables, motif_matrix),
    Selection_Frequency = rowMeans(!is.na(coefficient_matrix)),
    N_Selected_Runs = rowSums(!is.na(coefficient_matrix)),
    Mean_Coefficient = apply(coefficient_matrix, 1, safe_selected_mean),
    Median_Coefficient = apply(coefficient_matrix, 1, safe_selected_median),
    Combined_P_value = combined_p_values,
    Adjusted_P_value = adjusted_p_values,
    stringsAsFactors = FALSE
  )
  score_summary <- score_summary[
    order(score_summary$Adjusted_P_value, score_summary$Combined_P_value, na.last = TRUE),
  ]

  # DIAGNOSTIC: Check what variables were selected
  all_vars <- names(combined_p_values)
  motif_tfs <- rownames(motif_matrix)
  motif_cres <- colnames(motif_matrix)
  potential_tfs <- intersect(all_vars, motif_tfs)
  potential_cres <- intersect(all_vars, motif_cres)
  sig_tfs_vars <- potential_tfs[adjusted_p_values[potential_tfs] < p_value_cutoff & !is.na(adjusted_p_values[potential_tfs])]
  sig_cres_vars <- potential_cres[adjusted_p_values[potential_cres] < p_value_cutoff & !is.na(adjusted_p_values[potential_cres])]

  cat("\n  DIAGNOSTIC: Variable selection summary\n")
  cat("    Total variables with p-values:", length(all_vars), "\n")
  cat("    TFs in p-values:", length(potential_tfs), "\n")
  cat("    CREs in p-values:", length(potential_cres), "\n")
  cat("    Significant TFs (adjusted p <", p_value_cutoff, "):", length(sig_tfs_vars), "\n")
  cat("    Significant CREs (adjusted p <", p_value_cutoff, "):", length(sig_cres_vars), "\n")
  if (length(sig_cres_vars) == 0 && length(sig_tfs_vars) > 0) {
    cat("    WARNING: TFs are significant but NO CREs are significant!\n")
    cat("    This means no TFs can be identified (require both TF and CRE to be significant)\n")
  }

  # Find significant TFs using adjusted p-values
  sig_tfs <- find_significant_tfs(adjusted_p_values, motif_matrix, p_value_cutoff)

  # Get adjusted p-values for significant TFs
  sig_tf_pvals <- numeric(length(sig_tfs))
  names(sig_tf_pvals) <- sig_tfs

  for (i in seq_along(sig_tfs)) {
    tf <- sig_tfs[i]
    if (tf %in% names(adjusted_p_values)) {
      sig_tf_pvals[i] <- adjusted_p_values[tf]
    } else {
      sig_tf_pvals[i] <- NA_real_  # TF identified via CRE binding, no direct p-value
    }
  }

  # Filter: only keep TFs that are truly significant (adjusted p < cutoff) or have NA (via CRE binding)
  keep_idx <- is.na(sig_tf_pvals) | sig_tf_pvals < p_value_cutoff
  sig_tfs <- sig_tfs[keep_idx]
  sig_tf_pvals <- sig_tf_pvals[keep_idx]

  cat("✓ Found", length(sig_tfs), "significant TFs (BH-adjusted p-value <", p_value_cutoff, ")\n")

  # Count statistics
  n_sig_tfs <- length(sig_tfs_vars)  # TFs with significant p-values
  n_sig_cres <- length(sig_cres_vars)  # CREs with significant p-values
  n_identified_tfs <- length(sig_tfs)  # TFs identified (both TF and binding CRE are significant)

  # Get p-values for significant TFs and CREs
  sig_tf_pvals_list <- adjusted_p_values[sig_tfs_vars]
  sig_cre_pvals_list <- adjusted_p_values[sig_cres_vars]

  return(list(
    significant_tfs = sig_tfs,  # Identified TFs (both TF and CRE significant)
    significant_tf_pvals = sig_tf_pvals,
    significant_tfs_only = sig_tfs_vars,  # TFs with significant p-values
    significant_tfs_only_pvals = sig_tf_pvals_list,
    significant_cres = sig_cres_vars,  # CREs with significant p-values
    significant_cres_pvals = sig_cre_pvals_list,
    combined_p_values = combined_p_values,
    adjusted_p_values = adjusted_p_values,
    p_val_matrix = p_val_matrix,
    coefficient_matrix = coefficient_matrix,
    bootstrap_details = bootstrap_details,
    score_summary = score_summary,
    bootstrap_diagnostics = data.frame(
      run = seq_len(n_runs),
      elapsed_seconds = run_seconds,
      n_selected_variables = run_selected_vars
    ),
    n_sig_tfs = n_sig_tfs,
    n_sig_cres = n_sig_cres,
    n_identified_tfs = n_identified_tfs
  ))
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Load matrices from CSV files (period-specific)
#' @param input_dir Input directory containing matrices
#' @param target_gene Target gene symbol
#' @param period Period number
#' @param apply_batch_correction Whether batch correction was applied (affects subdirectory)
#' @return List with X, Y, Z, motif matrices, and batch_info
load_matrices <- function(input_dir, target_gene, period, apply_batch_correction = TRUE) {
  cat("Loading matrices from:", input_dir, "\n")

  # Preferred layout: period_x/[celltype_y/]target_gene
  period_dir <- file.path(input_dir, paste0("period_", period), target_gene)

  if (!dir.exists(period_dir)) {
    # Try cell-type directory layout: period_x/celltype_*/target_gene
    period_root <- file.path(input_dir, paste0("period_", period))
    if (dir.exists(period_root)) {
      celltype_dirs <- list.dirs(period_root, recursive = FALSE, full.names = TRUE)
      candidate_dirs <- file.path(celltype_dirs, target_gene)
      existing_candidates <- candidate_dirs[dir.exists(candidate_dirs)]
      if (length(existing_candidates) > 0) {
        period_dir <- existing_candidates[1]
        cat("  Using cell-type directory layout:", period_dir, "\n")
      }
    }
  }

  if (!dir.exists(period_dir)) {
    # Backward compatibility: legacy layout with correction suffix.
    correction_suffix <- if (apply_batch_correction) "batch_corrected" else "no_batch_correction"
    legacy_dir <- file.path(input_dir, paste0("period_", period), correction_suffix, target_gene)
    if (dir.exists(legacy_dir)) {
      period_dir <- legacy_dir
      cat("  Warning: Using legacy data-retrieval layout:", period_dir, "\n")
    } else {
      stop("ERROR: Input directory not found for period ", period, ", target ", target_gene)
    }
  }

  # Load X matrix (TF expression)
  X_file <- file.path(period_dir, paste0(target_gene, "_X_matrix.csv"))
  if (!file.exists(X_file)) {
    stop("ERROR: X matrix file not found: ", X_file)
  }
  X_matrix <- read.csv(X_file, row.names = 1, check.names = FALSE)
  cat("✓ Loaded X matrix:", dim(X_matrix), "\n")

  # Load Y matrix (target gene expression)
  Y_file <- file.path(period_dir, paste0(target_gene, "_Y_matrix.csv"))
  if (!file.exists(Y_file)) {
    stop("ERROR: Y matrix file not found: ", Y_file)
  }
  Y_matrix <- read.csv(Y_file, row.names = 1, check.names = FALSE)
  cat("✓ Loaded Y matrix:", dim(Y_matrix), "\n")

  # Load Z matrix (ATAC accessibility)
  Z_file <- file.path(period_dir, paste0(target_gene, "_Z_matrix.csv"))
  if (!file.exists(Z_file)) {
    stop("ERROR: Z matrix file not found: ", Z_file)
  }
  Z_matrix <- read.csv(Z_file, row.names = 1, check.names = FALSE)
  cat("✓ Loaded Z matrix:", dim(Z_matrix), "\n")

  # Prefer motif matrix colocated with retrieved X/Y/Z matrices
  motif_file <- file.path(period_dir, paste0(target_gene, "_motif_matrix.csv"))
  if (!file.exists(motif_file)) {
    motif_file <- file.path(PERIOD_MOTIF_RESULTS_DIR, paste0("period_", period), target_gene,
                            paste0(target_gene, "_motif_matrix.csv"))
  }
  if (!file.exists(motif_file)) {
    stop("ERROR: Motif matrix file not found: ", motif_file)
  }
  motif_matrix <- read.csv(motif_file, row.names = 1, check.names = FALSE)
  cat("✓ Loaded motif matrix:", dim(motif_matrix), "\n")

  # Try to load batch info
  batch_info_file <- file.path(period_dir, paste0(target_gene, "_batch_info.csv"))
  batch_info <- NULL
  if (file.exists(batch_info_file)) {
    batch_info <- read.csv(batch_info_file, stringsAsFactors = FALSE)
    cat("✓ Loaded batch information:", nrow(batch_info), "metacells\n")
  } else {
    cat("  Warning: Batch information file not found\n")
  }

  return(list(
    X_matrix = X_matrix,
    Y_matrix = Y_matrix,
    Z_matrix = Z_matrix,
    motif_matrix = motif_matrix,
    batch_info = batch_info
  ))
}

#' Verify matrix compatibility
#' @param matrices List of matrices
#' @return TRUE if compatible, stops with error if not
verify_matrices <- function(matrices) {
  cat("\nVerifying matrix compatibility...\n")

  X <- matrices$X_matrix
  Y <- matrices$Y_matrix
  Z <- matrices$Z_matrix
  motif <- matrices$motif_matrix

  # Check row alignment (metacells)
  x_metacells <- rownames(X)
  y_metacells <- rownames(Y)
  z_metacells <- rownames(Z)

  common_metacells <- intersect(intersect(x_metacells, y_metacells), z_metacells)

  if (length(common_metacells) == 0) {
    stop("ERROR: No common metacells found across X, Y, Z matrices")
  }

  cat("  Common metacells:", length(common_metacells), "\n")

  # Check TF name matching between X and motif matrix
  x_tfs <- colnames(X)
  motif_tfs <- rownames(motif)
  matched_tfs <- intersect(x_tfs, motif_tfs)

  cat("  TFs in X matrix:", length(x_tfs), "\n")
  cat("  TFs in motif matrix:", length(motif_tfs), "\n")
  cat("  Matched TFs:", length(matched_tfs), "\n")

  if (length(matched_tfs) == 0) {
    warning("WARNING: No TF names matched between X matrix and motif matrix!")
    cat("  Sample X matrix TFs:", paste(head(x_tfs, 5), collapse = ", "), "\n")
    cat("  Sample motif matrix TFs:", paste(head(motif_tfs, 5), collapse = ", "), "\n")
  }

  # Check peak name matching between Z and motif matrix
  z_peaks <- colnames(Z)
  motif_peaks <- colnames(motif)

  # Normalize peak names: convert dots to dashes for consistent matching
  # Z matrix might have dots (chr12.24469992.24470492) while motif has dashes (chr12-24469992-24470492)
  z_peaks_normalized <- gsub("\\.", "-", z_peaks)
  motif_peaks_normalized <- gsub("\\.", "-", motif_peaks)

  # Match using normalized names
  matched_peaks_normalized <- intersect(z_peaks_normalized, motif_peaks_normalized)

  # Get original peak names that match (use Z matrix names as reference)
  matched_peaks <- z_peaks[z_peaks_normalized %in% matched_peaks_normalized]

  cat("  Peaks in Z matrix:", length(z_peaks), "\n")
  cat("  Peaks in motif matrix:", length(motif_peaks), "\n")
  cat("  Matched peaks (after normalization):", length(matched_peaks), "\n")

  if (length(matched_peaks) == 0) {
    warning("WARNING: No peak names matched between Z matrix and motif matrix!")
    cat("  Sample Z matrix peaks:", paste(head(z_peaks, 3), collapse = ", "), "\n")
    cat("  Sample motif matrix peaks:", paste(head(motif_peaks, 3), collapse = ", "), "\n")
    cat("  Sample Z peaks (normalized):", paste(head(z_peaks_normalized, 3), collapse = ", "), "\n")
    cat("  Sample motif peaks (normalized):", paste(head(motif_peaks_normalized, 3), collapse = ", "), "\n")
  }

  # Subset matrices to common elements
  cat("\nSubsetting matrices to common elements...\n")

  # Subset to common metacells
  X_subset <- X[common_metacells, , drop = FALSE]
  Y_subset <- Y[common_metacells, , drop = FALSE]
  Z_subset <- Z[common_metacells, , drop = FALSE]

  # Subset batch info if available
  batch_info_subset <- NULL
  if (!is.null(matrices$batch_info)) {
    batch_info_subset <- matrices$batch_info[matrices$batch_info$metacell_id %in% common_metacells, , drop = FALSE]
    # Ensure order matches
    batch_info_subset <- batch_info_subset[match(common_metacells, batch_info_subset$metacell_id), , drop = FALSE]
  }

  # Subset motif matrix to matched TFs and peaks
  # Need to map matched_peaks (from Z) to motif matrix column names
  # matched_peaks are from Z matrix (might have dots), need to find corresponding motif peaks (with dashes)
  z_peaks_normalized <- gsub("\\.", "-", matched_peaks)
  motif_peak_indices <- match(z_peaks_normalized, colnames(motif))
  motif_peak_indices <- motif_peak_indices[!is.na(motif_peak_indices)]
  matched_motif_peaks <- colnames(motif)[motif_peak_indices]

  motif_subset <- motif[matched_tfs, matched_motif_peaks, drop = FALSE]

  # Subset X to matched TFs
  X_subset <- X_subset[, matched_tfs, drop = FALSE]

  # Subset Z to matched peaks (keep original Z peak names)
  Z_subset <- Z_subset[, matched_peaks, drop = FALSE]

  # Normalize Z matrix column names to use dashes (to match motif matrix)
  # This ensures consistency for downstream analysis
  colnames(Z_subset) <- gsub("\\.", "-", colnames(Z_subset))

  # Also normalize motif matrix column names to match
  colnames(motif_subset) <- gsub("\\.", "-", colnames(motif_subset))

  cat("  Final X matrix:", dim(X_subset), "\n")
  cat("  Final Y matrix:", dim(Y_subset), "\n")
  cat("  Final Z matrix:", dim(Z_subset), "\n")
  cat("  Final motif matrix:", dim(motif_subset), "\n")

  return(list(
    X_matrix = X_subset,
    Y_matrix = Y_subset,
    Z_matrix = Z_subset,
    motif_matrix = motif_subset,
    batch_info = batch_info_subset
  ))
}

#' Save inference results
#' @param results GRN inference results
#' @param target_gene Target gene symbol
#' @param period Period number
#' @param output_dir Output directory
#' @param p_value_cutoff P-value cutoff
#' @param apply_batch_correction Whether batch correction was applied
#' @param include_batch_covariates Whether batch covariates were included
save_results <- function(results, target_gene, period, output_dir, p_value_cutoff,
                         apply_batch_correction = TRUE, include_batch_covariates = FALSE) {
  # Flat output layout by period and target gene.
  period_output_dir <- file.path(output_dir, paste0("period_", period), target_gene)
  if (!dir.exists(period_output_dir)) {
    dir.create(period_output_dir, recursive = TRUE)
  }

  # Save three separate lists:
  # 1. Significant TFs (TF p-value < threshold)
  # 2. Significant CREs (CRE p-value < threshold)
  # 3. Identified TFs (both TF and binding CRE p-value < threshold)

  # 1. Save significant TFs (TF p-value < threshold)
  sig_tfs_only_file <- file.path(period_output_dir, paste0(target_gene, "_significant_tfs.csv"))
  if (length(results$significant_tfs_only) > 0) {
    sig_tfs_only_df <- data.frame(
      TF = results$significant_tfs_only,
      P_value = results$significant_tfs_only_pvals,
      stringsAsFactors = FALSE
    )
    sig_tfs_only_df <- sig_tfs_only_df[order(sig_tfs_only_df$P_value), ]
    write.csv(sig_tfs_only_df, sig_tfs_only_file, row.names = FALSE)
    cat("✓ Saved significant TFs to:", sig_tfs_only_file, "\n")
  } else {
    # Save blank CSV
    sig_tfs_only_df <- data.frame(
      TF = character(0),
      P_value = numeric(0),
      stringsAsFactors = FALSE
    )
    write.csv(sig_tfs_only_df, sig_tfs_only_file, row.names = FALSE)
    cat("✓ Saved blank significant TFs CSV to:", sig_tfs_only_file, "\n")
  }

  # 2. Save significant CREs (CRE p-value < threshold)
  sig_cres_file <- file.path(period_output_dir, paste0(target_gene, "_significant_cres.csv"))
  if (length(results$significant_cres) > 0) {
    sig_cres_df <- data.frame(
      CRE = results$significant_cres,
      P_value = results$significant_cres_pvals,
      stringsAsFactors = FALSE
    )
    sig_cres_df <- sig_cres_df[order(sig_cres_df$P_value), ]
    write.csv(sig_cres_df, sig_cres_file, row.names = FALSE)
    cat("✓ Saved significant CREs to:", sig_cres_file, "\n")
  } else {
    # Save blank CSV
    sig_cres_df <- data.frame(
      CRE = character(0),
      P_value = numeric(0),
      stringsAsFactors = FALSE
    )
    write.csv(sig_cres_df, sig_cres_file, row.names = FALSE)
    cat("✓ Saved blank significant CREs CSV to:", sig_cres_file, "\n")
  }

  # 3. Save identified TFs (both TF and binding CRE p-value < threshold)
  identified_tfs_file <- file.path(period_output_dir, paste0(target_gene, "_identified_tfs.csv"))
  sig_tf_file <- file.path(period_output_dir, paste0(target_gene, "_significant_tfs.txt"))

  if (length(results$significant_tfs) > 0) {
    # Create results dataframe and sort by p-value
    sig_tf_df <- data.frame(
      TF = results$significant_tfs,
      P_value = results$significant_tf_pvals,
      stringsAsFactors = FALSE
    )

    # Handle NA p-values (TFs identified via CRE binding but without direct p-values)
    sig_tf_df <- sig_tf_df[order(sig_tf_df$P_value, na.last = TRUE), ]

    # Write to text file with proper formatting
    writeLines(c(
      paste("Significant Transcription Factors for", target_gene),
      paste("Period:", period),
      paste(rep("=", 60), collapse = ""),
      paste("Target Gene:", target_gene),
      paste("Total Significant TFs:", length(results$significant_tfs)),
      paste("P-value Cutoff (BH-adjusted):", p_value_cutoff),
      paste(""),
      paste("Results:"),
      paste(rep("-", 60), collapse = ""),
      sprintf("%-25s %15s", "Transcription Factor", "P-value"),
      paste(rep("-", 60), collapse = "")
    ), sig_tf_file)

    # Append TF results
    for (i in seq_len(nrow(sig_tf_df))) {
      if (!is.na(sig_tf_df$P_value[i])) {
        write(sprintf("%-25s %15.4e", sig_tf_df$TF[i], sig_tf_df$P_value[i]),
              sig_tf_file, append = TRUE)
      } else {
        write(sprintf("%-25s %15s", sig_tf_df$TF[i], "NA (via CRE binding)"),
              sig_tf_file, append = TRUE)
      }
    }

    cat("✓ Saved significant TFs to:", sig_tf_file, "\n")

    # Save identified TFs as CSV
    write.csv(sig_tf_df, identified_tfs_file, row.names = FALSE)
    cat("✓ Saved identified TFs to:", identified_tfs_file, "\n")
  } else {
    # Write empty results file
    writeLines(c(
      paste("Identified Transcription Factors for", target_gene),
      paste("Period:", period),
      paste(rep("=", 60), collapse = ""),
      paste("Target Gene:", target_gene),
      paste("Total Identified TFs: 0"),
      paste("P-value Cutoff (BH-adjusted):", p_value_cutoff),
      paste(""),
      paste("No identified transcription factors found (requires both TF and binding CRE to be significant).")
    ), sig_tf_file)
    cat("✓ Saved results to:", sig_tf_file, "\n")

    # Save blank CSV file for identified TFs
    sig_tf_df <- data.frame(
      TF = character(0),
      P_value = numeric(0),
      stringsAsFactors = FALSE
    )
    write.csv(sig_tf_df, identified_tfs_file, row.names = FALSE)
    cat("✓ Saved blank identified TFs CSV to:", identified_tfs_file, "\n")
  }

  # Save combined p-values and adjusted p-values
  pval_df <- data.frame(
    Variable = names(results$combined_p_values),
    P_value = results$combined_p_values,
    Adjusted_P_value = results$adjusted_p_values,
    stringsAsFactors = FALSE
  )
  write.csv(pval_df, file.path(period_output_dir, paste0(target_gene, "_combined_p_values.csv")),
            row.names = FALSE)

  if (!is.null(results$bootstrap_diagnostics)) {
    diag_file <- file.path(period_output_dir, paste0(target_gene, "_bootstrap_diagnostics.csv"))
    write.csv(results$bootstrap_diagnostics, diag_file, row.names = FALSE)
    cat("✓ Saved bootstrap diagnostics to:", diag_file, "\n")
  }

  # Save summary
  summary_file <- file.path(period_output_dir, paste0(target_gene, "_inference_summary.txt"))

  # Get statistics (with defaults if not present for backward compatibility)
  n_sig_tfs <- ifelse(is.null(results$n_sig_tfs), length(results$significant_tfs), results$n_sig_tfs)
  n_sig_cres <- ifelse(is.null(results$n_sig_cres), 0, results$n_sig_cres)
  n_identified_tfs <- ifelse(is.null(results$n_identified_tfs), length(results$significant_tfs), results$n_identified_tfs)

  summary_lines <- c(
    paste("Period-Specific GRN Inference Pipeline Summary"),
    paste(rep("=", 60), collapse = ""),
    paste("Target Gene:", target_gene),
    paste("Period:", period),
    paste("Batch Correction (data retrieval):", ifelse(apply_batch_correction, "Yes", "No")),
    paste("Batch Covariates (in model):", ifelse(include_batch_covariates, "Yes", "No")),
    paste("Bootstrap Runs:", length(results$p_val_matrix[1,])),
    paste("P-value Cutoff (BH-adjusted):", p_value_cutoff),
    paste(""),
    paste("Results:"),
    paste("  Significant TFs (p-value < cutoff):", n_sig_tfs),
    paste("  Significant CREs (p-value < cutoff):", n_sig_cres),
    paste("  Identified TFs (TF and binding CRE both significant):", n_identified_tfs)
  )

  summary_lines <- c(summary_lines, "",
    paste("Output Files:"),
    paste("  -", paste0(target_gene, "_significant_tfs.csv"), "- Significant TFs (TF p-value < cutoff)"),
    paste("  -", paste0(target_gene, "_significant_cres.csv"), "- Significant CREs (CRE p-value < cutoff)"),
    paste("  -", paste0(target_gene, "_identified_tfs.csv"), "- Identified TFs (both TF and binding CRE p-value < cutoff)"))

  if (length(results$significant_tfs) == 0 && length(results$significant_tfs_only) == 0) {
    summary_lines <- c(summary_lines, "  No significant TFs found")
  }

  writeLines(summary_lines, summary_file)
}

# =============================================================================
# MAIN PIPELINE FUNCTION
# =============================================================================

#' Create batch dummy variables (one-hot encoding)
#' @param batch_labels Vector of batch labels for each sample
#' @return Matrix of batch dummy variables (samples × batches, one column per batch)
create_batch_dummies <- function(batch_labels) {
  unique_batches <- sort(unique(batch_labels))
  n_batches <- length(unique_batches)

  if (n_batches <= 1) {
    return(NULL)  # No need for batch dummies if only one batch
  }

  # Create one-hot encoding (but drop one level to avoid perfect collinearity)
  # Keep n_batches - 1 dummy variables
  batch_dummies <- matrix(0, nrow = length(batch_labels), ncol = n_batches - 1)
  colnames(batch_dummies) <- paste0("batch_", unique_batches[-1])

  for (i in 2:n_batches) {
    batch_dummies[, i - 1] <- as.numeric(batch_labels == unique_batches[i])
  }

  return(batch_dummies)
}

#' Run period-specific GRN inference pipeline
#' @param target_gene Target gene symbol
#' @param period Period number
#' @param n_runs Number of bootstrap runs (default: 100)
#' @param p_value_cutoff P-value cutoff (default: 0.05)
#' @param input_dir Input directory containing matrices (default: period_data_retrieval_results)
#' @param output_dir Output directory (default: period_grn_inference_results)
#' @param apply_batch_correction Whether batch correction was applied in data retrieval (affects input directory)
#' @param include_batch_covariates Whether to include batch indicators as unpenalized covariates (default: FALSE)
#' @return List with inference results
run_period_grn_inference_pipeline <- function(target_gene, period, n_runs = 100,
                                               p_value_cutoff = 0.05,
                                               input_dir = DEFAULT_INPUT_DIR,
                                               output_dir = DEFAULT_OUTPUT_DIR,
                                               apply_batch_correction = TRUE,
                                               include_batch_covariates = FALSE) {

  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("PERIOD-SPECIFIC GRN INFERENCE PIPELINE\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Target Gene:", target_gene, "\n")
  cat("Period:", period, "\n")
  cat("Bootstrap Runs:", n_runs, "\n")
  cat("P-value Cutoff:", p_value_cutoff, "\n")
  cat("Apply Batch Correction (data retrieval):", apply_batch_correction, "\n")
  cat("Include Batch Covariates:", include_batch_covariates, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")

  # Create output directory
  period_output_dir <- file.path(output_dir, paste0("period_", period), target_gene)
  if (!dir.exists(period_output_dir)) {
    dir.create(period_output_dir, recursive = TRUE)
  }

  # Step 1: Load matrices
  cat("Step 1: Loading matrices...\n")
  matrices <- load_matrices(input_dir, target_gene, period, apply_batch_correction)

  # Step 2: Verify and align matrices
  cat("\nStep 2: Verifying and aligning matrices...\n")
  matrices_aligned <- verify_matrices(matrices)

  # Safety guard: never include target-gene expression in predictor covariates (X).
  # Match case-insensitively so BCL11B/Bcl11B/bcl11b are treated as the same gene.
  x_cols_upper <- toupper(colnames(matrices_aligned$X_matrix))
  target_upper <- toupper(target_gene)
  target_in_x_idx <- which(x_cols_upper == target_upper)
  if (length(target_in_x_idx) > 0) {
    removed_cols <- colnames(matrices_aligned$X_matrix)[target_in_x_idx]
    matrices_aligned$X_matrix <- matrices_aligned$X_matrix[, -target_in_x_idx, drop = FALSE]
    cat("  Removed target gene from X covariates:", paste(removed_cols, collapse = ", "), "\n")
    cat("  X matrix after target-gene removal:", dim(matrices_aligned$X_matrix), "\n")
  } else {
    cat("  Target gene not present in X covariates (no removal needed)\n")
  }

  if (ncol(matrices_aligned$X_matrix) == 0) {
    stop("ERROR: X matrix has no predictors after removing target gene from covariates")
  }

  # Step 2.5: Create batch dummy variables if requested
  batch_covariates <- NULL
  if (include_batch_covariates) {
    if (!is.null(matrices_aligned$batch_info)) {
      # Get batch labels aligned with common metacells
      batch_labels <- matrices_aligned$batch_info$batch
      batch_covariates <- create_batch_dummies(batch_labels)

      if (!is.null(batch_covariates)) {
        cat("  Created", ncol(batch_covariates), "batch dummy variables\n")
      } else {
        cat("  Warning: Only one batch found, no batch covariates created\n")
      }
    } else {
      cat("  Warning: Batch information not available, cannot create batch covariates\n")
    }
  }

  # Step 3: Run GRN inference
  cat("\nStep 3: Running GRN inference...\n")
  progress_file <- file.path(period_output_dir, paste0(target_gene, "_bootstrap_progress.txt"))
  grn_results <- run_grn_inference(
    matrices_aligned$X_matrix,
    matrices_aligned$Y_matrix,
    matrices_aligned$Z_matrix,
    matrices_aligned$motif_matrix,
    batch_covariates = batch_covariates,
    n_runs = n_runs,
    p_value_cutoff = p_value_cutoff,
    progress_file = progress_file,
    progress_every = 5
  )

  # Step 4: Save results
  cat("\nStep 4: Saving results...\n")
  save_results(grn_results, target_gene, period, output_dir, p_value_cutoff,
               apply_batch_correction, include_batch_covariates)

  # Print final results
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("FINAL RESULTS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n")

  if (length(grn_results$significant_tfs) > 0) {
    cat("Significant TFs found:\n")
    for (i in seq_along(grn_results$significant_tfs)) {
      tf <- grn_results$significant_tfs[i]
      pval <- grn_results$significant_tf_pvals[i]
      if (!is.na(pval)) {
        cat("  -", tf, ":", sprintf("%.2e", pval), "\n")
      } else {
        cat("  -", tf, ": NA (via CRE binding)\n")
      }
    }
  } else {
    cat("No significant TFs found at p-value cutoff", p_value_cutoff, "\n")
  }

  cat("\nAnalysis completed successfully!\n")
  cat("Results saved to:", file.path(output_dir, paste0("period_", period), target_gene), "\n")

  return(grn_results)
}
