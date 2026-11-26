#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(grpregOverlap)
  library(progress)
})

classify_variable_type_locus <- function(vars, motif_matrix) {
  ifelse(
    vars %in% rownames(motif_matrix),
    "TF",
    ifelse(grepl("^period_", vars), "covariate", "CRE")
  )
}

safe_mean_na_locus <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}

safe_median_na_locus <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

run_grn_inference_with_scores_locus <- function(helpers,
                                                X_matrix,
                                                Y_matrix,
                                                Z_matrix,
                                                motif_matrix,
                                                batch_covariates = NULL,
                                                n_runs = 20,
                                                p_value_cutoff = 0.05,
                                                progress_file = NULL,
                                                progress_every = 5) {
  cat("\nRunning GRN inference with score capture...\n")

  TG_RNA <- as.numeric(Y_matrix[, 1])
  group <- helpers$construct_group(motif_matrix)
  cat("  Constructed", length(group), "groups for group lasso\n")

  if (length(group) == 0) {
    stop("No TF-CRE groups available after motif/data alignment")
  }

  if (!is.null(progress_file)) {
    writeLines(sprintf("0/%d", n_runs), progress_file)
  }

  pb <- progress::progress_bar$new(total = n_runs)
  bootstrap_rows <- vector("list", n_runs)
  run_seconds <- numeric(n_runs)
  run_selected_vars <- integer(n_runs)

  for (seed in seq_len(n_runs)) {
    run_t0 <- proc.time()[["elapsed"]]
    fit <- helpers$grpreg_fit(
      TF_RNA = X_matrix,
      CRE_ATAC = Z_matrix,
      TG_RNA = TG_RNA,
      group = group,
      batch_covariates = batch_covariates,
      verbose = 0,
      seed = seed
    )
    run_t1 <- proc.time()[["elapsed"]]
    run_seconds[seed] <- run_t1 - run_t0
    run_selected_vars[seed] <- ifelse(is.null(fit$n_selected_variables), NA_integer_, fit$n_selected_variables)

    coef_summary <- summary(fit$inference_model)$coefficients
    if (nrow(coef_summary) > 1) {
      raw_names <- rownames(coef_summary)[-1]
      normalized_names <- gsub("\\.", "-", raw_names)
      bootstrap_rows[[seed]] <- data.frame(
        run = seed,
        variable = normalized_names,
        variable_type = classify_variable_type_locus(normalized_names, motif_matrix),
        coefficient = coef_summary[-1, 1],
        post_selection_p_value = coef_summary[-1, 4],
        stringsAsFactors = FALSE
      )
    } else {
      bootstrap_rows[[seed]] <- data.frame(
        run = integer(0),
        variable = character(0),
        variable_type = character(0),
        coefficient = numeric(0),
        post_selection_p_value = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    pb$tick()
    if ((seed %% progress_every) == 0 || seed == 1 || seed == n_runs) {
      cat(sprintf("  Bootstrap progress: %d/%d\n", seed, n_runs))
    }
    cat(sprintf("    run=%d elapsed=%.3fs selected_vars=%d\n", seed, run_seconds[seed], run_selected_vars[seed]))
    if (!is.null(progress_file)) {
      writeLines(sprintf("%d/%d", seed, n_runs), progress_file)
    }
  }

  bootstrap_details <- do.call(rbind, bootstrap_rows)
  all_variables <- unique(bootstrap_details$variable)
  if (length(all_variables) == 0) {
    stop("No variables were selected in any bootstrap run")
  }

  p_val_matrix <- matrix(0.95, nrow = length(all_variables), ncol = n_runs)
  coef_matrix <- matrix(NA_real_, nrow = length(all_variables), ncol = n_runs)
  rownames(p_val_matrix) <- all_variables
  rownames(coef_matrix) <- all_variables

  for (seed in seq_len(n_runs)) {
    run_df <- bootstrap_rows[[seed]]
    if (nrow(run_df) == 0) next
    idx <- match(run_df$variable, all_variables)
    p_val_matrix[idx, seed] <- run_df$post_selection_p_value
    coef_matrix[idx, seed] <- run_df$coefficient
  }

  combined_p_values <- apply(p_val_matrix, 1, helpers$cauchy_combine_p)
  adjusted_p_values <- p.adjust(combined_p_values, method = "BH")

  score_summary <- data.frame(
    Variable = all_variables,
    Variable_Type = classify_variable_type_locus(all_variables, motif_matrix),
    Selection_Frequency = rowMeans(!is.na(coef_matrix)),
    N_Selected_Runs = rowSums(!is.na(coef_matrix)),
    Mean_Coefficient = apply(coef_matrix, 1, safe_mean_na_locus),
    Median_Coefficient = apply(coef_matrix, 1, safe_median_na_locus),
    Combined_P_value = combined_p_values,
    Adjusted_P_value = adjusted_p_values,
    stringsAsFactors = FALSE
  )
  score_summary <- score_summary[order(score_summary$Adjusted_P_value, score_summary$Combined_P_value, na.last = TRUE), ]

  tf_rows <- score_summary$Variable_Type == "TF"
  cre_rows <- score_summary$Variable_Type == "CRE"
  sig_tfs_only_df <- score_summary[tf_rows & !is.na(score_summary$Adjusted_P_value) & score_summary$Adjusted_P_value < p_value_cutoff, , drop = FALSE]
  sig_cres_df <- score_summary[cre_rows & !is.na(score_summary$Adjusted_P_value) & score_summary$Adjusted_P_value < p_value_cutoff, , drop = FALSE]

  sig_tfs <- helpers$find_significant_tfs(adjusted_p_values, motif_matrix, p_value_cutoff)
  identified_tfs_df <- score_summary[match(sig_tfs, score_summary$Variable), , drop = FALSE]
  identified_tfs_df <- identified_tfs_df[!is.na(identified_tfs_df$Variable), , drop = FALSE]

  cat(sprintf("  Bootstrap timing summary (seconds): min=%.3f median=%.3f mean=%.3f max=%.3f\n",
              min(run_seconds), median(run_seconds), mean(run_seconds), max(run_seconds)))
  cat(sprintf("  Selected-variable summary: min=%d median=%.1f mean=%.2f max=%d\n",
              min(run_selected_vars), median(run_selected_vars), mean(run_selected_vars), max(run_selected_vars)))

  list(
    bootstrap_details = bootstrap_details,
    bootstrap_diagnostics = data.frame(
      run = seq_len(n_runs),
      elapsed_seconds = run_seconds,
      n_selected_variables = run_selected_vars,
      stringsAsFactors = FALSE
    ),
    coefficient_matrix = coef_matrix,
    p_val_matrix = p_val_matrix,
    score_summary = score_summary,
    significant_tfs_only_df = sig_tfs_only_df,
    significant_cres_df = sig_cres_df,
    identified_tfs_df = identified_tfs_df,
    significant_tfs_only = sig_tfs_only_df$Variable,
    significant_tfs_only_pvals = sig_tfs_only_df$Adjusted_P_value,
    significant_cres = sig_cres_df$Variable,
    significant_cres_pvals = sig_cres_df$Adjusted_P_value,
    significant_tfs = identified_tfs_df$Variable,
    significant_tf_pvals = identified_tfs_df$Adjusted_P_value,
    combined_p_values = combined_p_values,
    adjusted_p_values = adjusted_p_values,
    n_sig_tfs = nrow(sig_tfs_only_df),
    n_sig_cres = nrow(sig_cres_df),
    n_identified_tfs = nrow(identified_tfs_df)
  )
}
