#!/usr/bin/env Rscript
# =====================================================================
#  coloc-SuSiE fine-mapping + colocalization pipeline
#  - Harmonizes GWAS, eQTL, and LD to a common reference allele
#  - Uses coloc::runsusie (canonical workflow) then coloc.susie
#  - Writes z-vs-LD diagnostics (estimate_s_rss) into the output
#  - Logs target loci (e.g. AHI1) with explicit drop reasons
#  - Exports a multi-sheet .xlsx of results
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(coloc)
  library(susieR)
  library(openxlsx)
})

## ============================ CONFIG ============================
# Cell type and stage are read from environment variables when set (by the
# launcher script), otherwise fall back to these defaults for standalone runs.
CELLTYPE <- Sys.getenv("CT",    unset = "ExN")        # ExN | InN | Glia | NNC
STAGE    <- Sys.getenv("STAGE", unset = "prenatal")   # prenatal | postnatal
TRAIT    <- Sys.getenv("TRAIT", unset = "MS")         # GWAS trait prefix

stopifnot(CELLTYPE %in% c("ExN","InN","Glia","NNC"),
          STAGE    %in% c("prenatal","postnatal"))

BASE_INPUT  <- "/mnt/s3/FXN/QTL/eQTL/coloc_susie"
INPUT_DIR   <- file.path(BASE_INPUT, paste0(TRAIT, "_bmind_", STAGE, "_", CELLTYPE))
LD_DIR      <- "/mnt/s3/FXN/wgs_bv1_and_bv2/upload/LD_exact_blocks"
N_GWAS      <- 41505
# eQTL sample size is stage-level (same across cell types):
#   prenatal = 127, postnatal = 100
N_EQTL_BY_STAGE <- c(prenatal = 127, postnatal = 100)
N_EQTL <- N_EQTL_BY_STAGE[[STAGE]]
if (is.null(N_EQTL) || is.na(N_EQTL))
  stop("N_EQTL not set for stage '", STAGE, "'")
TAG         <- paste0(TRAIT, "_", CELLTYPE, "_", STAGE)   # output name stem
OUT_TSV     <- paste0(TAG, "_coloc_results.tsv")
OUT_XLSX    <- paste0(TAG, "_coloc_results.xlsx")
OUT_RDS     <- paste0(TAG, "_coloc_results.rds")
MAX_CS      <- 5
MIN_OVERLAP <- 5        # minimum SNPs shared across GWAS/eQTL/LD
L_SUSIE     <- 5        # max causal signals per locus (10 over-splits at small N)
MIN_ABS_CORR<- 0.5      # CS purity filter: drop credible sets with min|r| below this
                        # (0.5 is the susieR default; raise toward 0.8 to be stricter)
COVERAGE    <- 0.95
DEBUG       <- TRUE

# Loci you KNOW should appear; the pipeline logs exactly why if they drop.
TARGET_LOCI <- c("AHI1")

# When TRUE, run only the diagnostic block (file counts, column names,
# LD coverage, target-locus check) and STOP before the main pipeline.
# Set to FALSE once diagnostics look right, then re-run for the full analysis.
DIAGNOSE_ONLY <- as.logical(Sys.getenv("DIAGNOSE_ONLY", unset = "FALSE"))
## ================================================================

message("=== RUN CONFIG: ", TAG, " ===")
message("    INPUT_DIR = ", INPUT_DIR)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# ---------------------------------------------------------------
#  LD block index
# ---------------------------------------------------------------
load_ld_blocks <- function(ld_dir) {
  files <- list.files(ld_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) return(tibble(file = character(0), chr = character(0)))
  map_df(files, function(f) {
    base <- gsub("_hg19\\.rds$|\\.rds$", "", basename(f))
    tibble(file = f, chr = strsplit(base, "[:_]")[[1]][1])
  })
}

# ---------------------------------------------------------------
#  Allele harmonization
#  Aligns GWAS beta, eQTL beta, and the LD matrix to a single
#  reference orientation per SNP. Returns aligned vectors + LD.
#  This is the step that prevents allele-swapped loci (e.g. AHI1)
#  from being dropped or sign-flipped relative to LD.
# ---------------------------------------------------------------
harmonize <- function(df, ld, debug = TRUE) {
  # df carries effect/non-effect alleles per source:
  #   eff_gwas, noneff_gwas, beta_gwas
  #   eff_eqtl, noneff_eqtl, beta_eqtl
  # Goal: express the eQTL effect on the SAME allele as the GWAS effect,
  # so both betas (and the LD R, which encodes the GWAS reference) agree.
  flip <- function(a) chartr("ACGT", "TGCA", a)

  # eQTL orientation relative to GWAS effect allele
  same  <- df$eff_eqtl == df$eff_gwas & df$noneff_eqtl == df$noneff_gwas
  swap  <- df$eff_eqtl == df$noneff_gwas & df$noneff_eqtl == df$eff_gwas
  fsame <- flip(df$eff_eqtl) == df$eff_gwas & flip(df$noneff_eqtl) == df$noneff_gwas
  fswap <- flip(df$eff_eqtl) == df$noneff_gwas & flip(df$noneff_eqtl) == df$eff_gwas

  resolvable <- same | swap | fsame | fswap
  # When eQTL effect is on the opposite allele, flip its beta to GWAS orientation
  df$beta_eqtl[swap | fswap] <- -df$beta_eqtl[swap | fswap]

  # Drop strand-ambiguous A/T and C/G SNPs (cannot resolve orientation safely)
  ambiguous <- (df$eff_gwas == flip(df$noneff_gwas))
  keep <- resolvable & !ambiguous &
          !is.na(df$beta_gwas) & !is.na(df$beta_eqtl) &
          !is.na(df$se_gwas) & df$se_gwas > 0 &
          !is.na(df$se_eqtl) & df$se_eqtl > 0

  if (debug) {
    message("    harmonize: ", sum(keep), "/", nrow(df),
            " SNPs kept (", sum((swap | fswap) & keep), " eQTL-beta-flipped, ",
            sum(ambiguous), " strand-ambiguous dropped, ",
            sum(!resolvable), " allele-mismatch dropped)")
  }
  df[keep, , drop = FALSE]
}

# ---------------------------------------------------------------
#  Credible-set extraction from a runsusie object
#  runsusie stores snp names, so $sets$cs entries are named/indexed
#  consistently; we read coverage & purity directly from $sets.
# ---------------------------------------------------------------
extract_cs <- function(s, prefix, max_cs) {
  out <- list()
  out[[paste0(prefix, "_n_cs")]] <- 0L
  for (i in seq_len(max_cs)) {
    out[[paste0(prefix, "_cs", i, "_snps")]]     <- NA_character_
    out[[paste0(prefix, "_cs", i, "_size")]]     <- 0L
    out[[paste0(prefix, "_cs", i, "_top_snp")]]  <- NA_character_
    out[[paste0(prefix, "_cs", i, "_top_pip")]]  <- NA_real_
    out[[paste0(prefix, "_cs", i, "_purity")]]   <- NA_real_
    out[[paste0(prefix, "_cs", i, "_coverage")]] <- NA_real_
  }
  if (is.null(s) || is.null(s$sets) || is.null(s$sets$cs)) return(out)

  cs   <- s$sets$cs
  pip  <- s$pip
  pur  <- s$sets$purity
  cov  <- s$sets$coverage
  out[[paste0(prefix, "_n_cs")]] <- length(cs)

  for (i in seq_along(cs)) {
    if (i > max_cs) break
    idx  <- cs[[i]]
    snps <- names(pip)[idx]
    p    <- pip[idx]
    ord  <- order(-p)
    out[[paste0(prefix, "_cs", i, "_snps")]]    <- paste(snps[ord], collapse = ",")
    out[[paste0(prefix, "_cs", i, "_size")]]    <- length(idx)
    out[[paste0(prefix, "_cs", i, "_top_snp")]] <- snps[ord][1]
    out[[paste0(prefix, "_cs", i, "_top_pip")]] <- round(p[ord][1], 4)
    if (!is.null(pur) && nrow(pur) >= i)
      out[[paste0(prefix, "_cs", i, "_purity")]] <- round(pur[i, "min.abs.corr"], 4)
    if (!is.null(cov) && length(cov) >= i)
      out[[paste0(prefix, "_cs", i, "_coverage")]] <- round(cov[i], 4)
  }
  out
}

# ---------------------------------------------------------------
#  Per-locus coloc-SuSiE
# ---------------------------------------------------------------
run_one <- function(ob, ld, chr_name, gene_ensembl, gene_name,
                    N_gwas, N_eqtl, max_cs, debug = TRUE) {

  base <- list(
    chr = chr_name, gene_name = gene_name, gene_ensembl = gene_ensembl,
    nsnps = NA_integer_, status = NA_character_, method = NA_character_,
    PP_H0 = NA_real_, PP_H1 = NA_real_, PP_H2 = NA_real_,
    PP_H3 = NA_real_, PP_H4 = NA_real_,
    best_H4_snp = NA_character_, best_H4 = NA_real_,
    gwas_lambda = NA_real_, eqtl_lambda = NA_real_   # z-vs-LD diagnostic
  )
  fail <- function(reason) { base$status <- reason; as.data.frame(base) }

  if (nrow(ob) < MIN_OVERLAP) return(fail("too_few_overlap"))

  # ---- clean inputs so runsusie/susie_rss never sees NA/Inf ----
  # 1. require finite beta & positive finite SE in both datasets
  ok <- is.finite(ob$beta_gwas) & is.finite(ob$se_gwas) & ob$se_gwas > 0 &
        is.finite(ob$beta_eqtl) & is.finite(ob$se_eqtl) & ob$se_eqtl > 0 &
        is.finite(ob$MAF) & ob$MAF > 0 & ob$MAF < 1
  ob <- ob[ok, , drop = FALSE]
  ld <- ld[ob$SNP, ob$SNP, drop = FALSE]

  # 2. drop SNPs participating in any NA/non-finite LD cell, iteratively
  if (any(!is.finite(ld))) {
    repeat {
      bad_count <- rowSums(!is.finite(ld))
      if (all(bad_count == 0)) break
      drop1 <- which.max(bad_count)        # remove worst offender, re-check
      ld <- ld[-drop1, -drop1, drop = FALSE]
      ob <- ob[-drop1, , drop = FALSE]
      if (nrow(ld) < MIN_OVERLAP) break
    }
  }
  # 3. force exact symmetry and unit diagonal (guards tiny asymmetries)
  ld <- (ld + t(ld)) / 2
  diag(ld) <- 1

  # 4. PSD-repair once, ourselves, so susieR doesn't do it repeatedly and
  #    leave a subtly inconsistent matrix. Clip negative eigenvalues to a
  #    small positive floor, rescale back to a correlation matrix.
  ev <- tryCatch(eigen(ld, symmetric = TRUE), error = function(e) NULL)
  if (!is.null(ev) && min(ev$values) < 1e-4) {
    ev$values[ev$values < 1e-4] <- 1e-4
    ld <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
    d  <- sqrt(diag(ld)); ld <- ld / outer(d, d)   # back to correlation scale
    ld <- (ld + t(ld)) / 2; diag(ld) <- 1
    dimnames(ld) <- list(ob$SNP, ob$SNP)
  }
  if (any(!is.finite(ld))) return(fail("LD_not_repairable"))

  if (nrow(ob) < MIN_OVERLAP) return(fail("too_few_after_LD_clean"))

  # z-scores aligned to the (already harmonized & ordered) LD
  z_gwas <- ob$beta_gwas / ob$se_gwas
  z_eqtl <- ob$beta_eqtl / ob$se_eqtl
  base$nsnps <- nrow(ob)

  # z-vs-LD consistency diagnostic (flags misaligned LD/alleles)
  base$gwas_lambda <- tryCatch(
    estimate_s_rss(z = z_gwas, R = ld, n = N_gwas), error = function(e) NA_real_)
  base$eqtl_lambda <- tryCatch(
    estimate_s_rss(z = z_eqtl, R = ld, n = N_eqtl), error = function(e) NA_real_)

  # ---- detect & correct per-SNP allele flips vs the LD panel ----
  # kriging_rss predicts each z from the others via LD; SNPs whose observed
  # and predicted z have strongly opposite sign are flipped relative to LD.
  # We flip beta (and z) so the summary stats match the LD orientation.
  fix_flips <- function(z, beta, n, tag) {
    kr <- tryCatch(kriging_rss(z = z, R = ld, n = n), error = function(e) NULL)
    if (is.null(kr) || is.null(kr$conditional_dist)) return(list(beta = beta, z = z, nflip = 0))
    cd <- kr$conditional_dist
    logLR <- cd$logLR %||% rep(0, length(z))
    pred  <- cd$conditional_mean %||% rep(0, length(z))
    # NA-proof: treat any non-finite diagnostic as "no flip" for that SNP
    logLR[!is.finite(logLR)] <- 0
    pred[!is.finite(pred)]   <- 0
    cond <- logLR > 2 & sign(pred) == -sign(z) & abs(z) > 2
    cond[is.na(cond)] <- FALSE
    flip <- which(cond)
    if (length(flip) > 0) {
      beta[flip] <- -beta[flip]; z[flip] <- -z[flip]
      if (debug) message("    ", tag, ": flipped ", length(flip),
                         " SNP(s) to match LD orientation")
    }
    list(beta = beta, z = z, nflip = length(flip))
  }
  fg <- fix_flips(z_gwas, ob$beta_gwas, N_gwas, "GWAS")
  fe <- fix_flips(z_eqtl, ob$beta_eqtl, N_eqtl, "eQTL")
  ob$beta_gwas <- fg$beta; ob$beta_eqtl <- fe$beta

  # coloc dataset objects
  D_gwas <- list(beta = ob$beta_gwas, varbeta = ob$se_gwas^2,
                 snp = ob$SNP, position = ob$BP, type = "quant",
                 N = N_gwas, sdY = 1, LD = ld)
  D_eqtl <- list(beta = ob$beta_eqtl, varbeta = ob$se_eqtl^2,
                 snp = ob$SNP, position = ob$BP, type = "quant",
                 N = N_eqtl, MAF = ob$MAF, LD = ld)

  if (!is.null(check_dataset(D_gwas, req = "LD"))) return(fail("bad_gwas_dataset"))
  if (!is.null(check_dataset(D_eqtl, req = "LD"))) return(fail("bad_eqtl_dataset"))

  # Your LD is an EXTERNAL panel (WGS), not in-sample eQTL/GWAS genotypes.
  # estimate_residual_variance = TRUE assumes in-sample LD and fails
  # ("estimated value is negative") with external LD, so we keep it FALSE.
  # Overfitting (spurious size-1 credible sets) is instead controlled by:
  #   - a modest L (set in CONFIG)
  #   - SuSiE's built-in purity filter: credible sets with min |r| below
  #     min_abs_corr are dropped, which removes the spurious singletons.
  s_gwas <- tryCatch(
    runsusie(D_gwas, L = L_SUSIE, coverage = COVERAGE, maxit = 1000,
             min_abs_corr = MIN_ABS_CORR,
             check_prior = FALSE, estimate_residual_variance = FALSE,
             repeat_until_convergence = TRUE),
    error = function(e) { if (debug) message("    runsusie GWAS: ", e$message); NULL })
  s_eqtl <- tryCatch(
    runsusie(D_eqtl, L = L_SUSIE, coverage = COVERAGE, maxit = 1000,
             min_abs_corr = MIN_ABS_CORR,
             check_prior = FALSE, estimate_residual_variance = FALSE,
             repeat_until_convergence = TRUE),
    error = function(e) { if (debug) message("    runsusie eQTL: ", e$message); NULL })

  res <- list(c(base, extract_cs(s_gwas, "gwas", max_cs),
                extract_cs(s_eqtl, "eqtl", max_cs)))[[1]]

  # coloc.abf fallback: single-causal-variant, LD-free, always returns a
  # posterior. Used when coloc.susie yields no credible-set pair.
  run_abf <- function() {
    da <- list(beta = ob$beta_gwas, varbeta = ob$se_gwas^2, snp = ob$SNP,
               position = ob$BP, type = "quant", N = N_gwas, sdY = 1)
    db <- list(beta = ob$beta_eqtl, varbeta = ob$se_eqtl^2, snp = ob$SNP,
               position = ob$BP, type = "quant", N = N_eqtl, MAF = ob$MAF)
    cr <- tryCatch(coloc.abf(da, db),
                   error = function(e) { if (debug) message("    coloc.abf: ", e$message); NULL })
    if (is.null(cr) || is.null(cr$summary)) return(FALSE)
    s <- cr$summary
    res$PP_H0 <<- round(as.numeric(s[["PP.H0.abf"]]), 4)
    res$PP_H1 <<- round(as.numeric(s[["PP.H1.abf"]]), 4)
    res$PP_H2 <<- round(as.numeric(s[["PP.H2.abf"]]), 4)
    res$PP_H3 <<- round(as.numeric(s[["PP.H3.abf"]]), 4)
    res$PP_H4 <<- round(as.numeric(s[["PP.H4.abf"]]), 4)
    if (!is.null(cr$results) && "SNP.PP.H4" %in% names(cr$results)) {
      rd <- cr$results; bi <- which.max(rd$SNP.PP.H4)
      res$best_H4_snp <<- as.character(rd$snp[bi] %||% NA_character_)
      res$best_H4     <<- round(rd$SNP.PP.H4[bi], 4)
    }
    TRUE
  }

  pairs_df <- NULL   # long-format: one row per colocalizing CS-pair

  if (!is.null(s_gwas) && !is.null(s_eqtl)) {
    cr <- tryCatch(coloc.susie(s_gwas, s_eqtl),
                   error = function(e) { if (debug) message("    coloc.susie: ", e$message); NULL })
    if (!is.null(cr) && !is.null(cr$summary) && nrow(cr$summary) > 0) {
      sm <- as.data.frame(cr$summary)
      best <- which.max(sm$PP.H4.abf)
      res$PP_H0 <- round(sm$PP.H0.abf[best], 4)
      res$PP_H1 <- round(sm$PP.H1.abf[best], 4)
      res$PP_H2 <- round(sm$PP.H2.abf[best], 4)
      res$PP_H3 <- round(sm$PP.H3.abf[best], 4)
      res$PP_H4 <- round(sm$PP.H4.abf[best], 4)

      # ---- build long-format pair table (standard supplement form) ----
      # coloc's idx1/idx2 are credible-set LABELS (matching names of $sets$cs,
      # e.g. "L1","L3"), not positional indices. Resolve by name, then position,
      # with bounds checks so a label mismatch never crashes the locus.
      lead_of <- function(s, idx) {
        cs_list <- s$sets$cs
        if (is.null(cs_list) || length(cs_list) == 0) return(c(snp = NA, pip = NA))
        cs <- NULL
        nm <- names(cs_list)
        # try name match: "L3" or "3"
        cand <- c(as.character(idx), paste0("L", idx))
        hit  <- which(nm %in% cand)
        if (length(hit) >= 1) cs <- cs_list[[hit[1]]]
        else if (is.numeric(idx) && idx >= 1 && idx <= length(cs_list)) cs <- cs_list[[idx]]
        if (is.null(cs) || length(cs) == 0) return(c(snp = NA, pip = NA))
        cs <- cs[cs >= 1 & cs <= length(s$pip)]
        if (length(cs) == 0) return(c(snp = NA, pip = NA))
        pip <- s$pip[cs]; k <- which.max(pip)
        c(snp = names(s$pip)[cs][k], pip = round(unname(pip[k]), 4))
      }
      gi <- if ("idx1" %in% names(sm)) sm$idx1 else seq_len(nrow(sm))
      ei <- if ("idx2" %in% names(sm)) sm$idx2 else seq_len(nrow(sm))
      pairs_df <- tryCatch(data.frame(
        gene_name    = gene_name, chr = chr_name,
        gwas_cs      = gi, eqtl_cs = ei,
        gwas_lead    = sapply(gi, function(k) lead_of(s_gwas, k)["snp"]),
        gwas_lead_pip= sapply(gi, function(k) lead_of(s_gwas, k)["pip"]),
        eqtl_lead    = sapply(ei, function(k) lead_of(s_eqtl, k)["snp"]),
        eqtl_lead_pip= sapply(ei, function(k) lead_of(s_eqtl, k)["pip"]),
        PP_H0 = round(sm$PP.H0.abf, 4), PP_H1 = round(sm$PP.H1.abf, 4),
        PP_H2 = round(sm$PP.H2.abf, 4), PP_H3 = round(sm$PP.H3.abf, 4),
        PP_H4 = round(sm$PP.H4.abf, 4),
        stringsAsFactors = FALSE, row.names = NULL),
        error = function(e) { if (debug) message("    pairs table skipped: ", e$message); NULL })

      if (!is.null(cr$results)) {
        rd <- as.data.frame(cr$results)
        h4 <- grep("SNP.PP.H4", names(rd), value = TRUE)
        if (length(h4) >= 1) {
          bi <- which.max(rd[[h4[1]]])
          res$best_H4_snp <- as.character(rd$snp[bi] %||% NA_character_)
          res$best_H4     <- round(rd[[h4[1]]][bi], 4)
        }
      }
      res$status <- "ok"
      res$method <- "susie"
    } else {
      if (run_abf()) {
        res$status <- "ok"
        res$method <- "abf_fallback_no_cs"
      } else {
        res$status <- "coloc_no_cs"
        res$method <- "none"
      }
    }
  } else {
    if (run_abf()) {
      res$status <- "ok"
      res$method <- "abf_fallback_susie_failed"
    } else {
      res$status <- "susie_failed"
      res$method <- "none"
    }
  }
  out_row <- as.data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  attr(out_row, "pairs") <- pairs_df
  out_row
}

# ---------------------------------------------------------------
#  Pipeline
# ---------------------------------------------------------------
run_pipeline <- function() {
  ld_blocks  <- load_ld_blocks(LD_DIR)
  gwas_files <- sort(list.files(INPUT_DIR, pattern = "gwas", full.names = TRUE))
  eqtl_files <- sort(list.files(INPUT_DIR, pattern = "eqtl", full.names = TRUE))
  stopifnot(length(gwas_files) > 0, length(eqtl_files) > 0)
  message("Found ", length(gwas_files), " GWAS and ", length(eqtl_files), " eQTL files")

  all_results <- list()
  all_pairs   <- list()
  target_log  <- list()

  for (i in seq_along(gwas_files)) {
    gwas <- suppressMessages(readr::read_tsv(gwas_files[i], show_col_types = FALSE)) |> as.data.frame()
    eqtl <- suppressMessages(readr::read_tsv(eqtl_files[min(i, length(eqtl_files))], show_col_types = FALSE)) |> as.data.frame()

    gene_name <- if ("GENESYMBOL" %in% names(eqtl)) unique(eqtl$GENESYMBOL)[1] else paste0("gene_", i)
    gene_ens  <- if ("EnsemblID"  %in% names(eqtl)) unique(eqtl$EnsemblID)[1]  else paste0("unknown_", i)
    is_target <- gene_name %in% TARGET_LOCI
    log_target <- function(reason, n = NA) if (is_target)
      target_log[[gene_name]] <<- tibble(gene = gene_name, reason = reason, n_snps = n)

    message("\n=== [", i, "/", length(gwas_files), "] ", gene_name,
            if (is_target) "  <-- TARGET LOCUS" else "")

    # normalize key columns to character / numeric
    for (col in c("CHR","A1","A2","SNP")) {
      if (col %in% names(gwas)) gwas[[col]] <- toupper(as.character(gwas[[col]]))
      if (col %in% names(eqtl)) eqtl[[col]] <- toupper(as.character(eqtl[[col]]))
    }
    # ---- resolve column names flexibly (adapts to your file headers) ----
    pick <- function(df, candidates, what, required = TRUE) {
      hit <- candidates[candidates %in% names(df)]
      if (length(hit) == 0) {
        if (required)
          stop("Could not find a '", what, "' column. Looked for: ",
               paste(candidates, collapse = "/"),
               " | available: ", paste(names(df), collapse = ", "))
        return(NA_character_)
      }
      hit[1]
    }
    g_snp  <- pick(gwas, c("SNP","rsid","RSID","variant_id","ID"), "GWAS SNP")
    g_a1   <- pick(gwas, c("A1","ALT","EA","effect_allele","Allele1"), "GWAS A1")
    g_a2   <- pick(gwas, c("A2","REF","NEA","other_allele","Allele2"), "GWAS A2")
    g_se   <- pick(gwas, c("SE","se","StdErr","standard_error","sebeta"), "GWAS SE", required = FALSE)
    g_beta <- pick(gwas, c("BETA","beta","b","Effect","effect"), "GWAS beta", required = FALSE)
    g_or   <- pick(gwas, c("OR","odds_ratio"), "GWAS OR", required = FALSE)
    g_p    <- pick(gwas, c("P","p","pval","P_value","pvalue"), "GWAS P", required = FALSE)
    g_frq  <- pick(gwas, c("FRQ","FREQ","MAF","EAF","A1FREQ","af"), "GWAS FRQ", required = FALSE)
    g_bp   <- pick(gwas, c("BP","POS","position","pos","BPos"), "GWAS BP")
    g_chr  <- pick(gwas, c("CHR","chr","chromosome","Chr"), "GWAS CHR")

    e_snp  <- pick(eqtl, c("SNP","rsid","RSID","variant_id","ID"), "eQTL SNP")
    e_a1   <- pick(eqtl, c("A1","ALT","EA","effect_allele","Allele1"), "eQTL A1")
    e_a2   <- pick(eqtl, c("A2","REF","NEA","other_allele","Allele2"), "eQTL A2")
    e_eff  <- pick(eqtl, c("alt_alleles","ALT","effect_allele","EA"), "eQTL effect allele", required = FALSE)
    e_se   <- pick(eqtl, c("SE","se","StdErr","standard_error","sebeta"), "eQTL SE")
    e_beta <- pick(eqtl, c("BETA","beta","b","slope","Effect"), "eQTL beta")
    e_frq  <- pick(eqtl, c("FRQ","FREQ","MAF","EAF","af","AF"), "eQTL FRQ", required = FALSE)

    chr_name <- paste0("chr", gsub("^chr", "", unique(gwas[[g_chr]])[1], ignore.case = TRUE))

    # ---- GWAS beta/SE: reconstruct SE from OR + P when SE is absent ----
    # beta = log(OR); z = sign(beta)*qnorm(P/2, lower=FALSE); SE = beta/z
    if (!is.na(g_beta)) {
      gb <- safe_num(gwas[[g_beta]])
    } else if (!is.na(g_or)) {
      gb <- log(safe_num(gwas[[g_or]]))
    } else {
      stop("GWAS has neither BETA nor OR")
    }
    if (!is.na(g_se)) {
      gse <- safe_num(gwas[[g_se]])
    } else if (!is.na(g_p)) {
      gp  <- safe_num(gwas[[g_p]])
      gp  <- pmin(pmax(gp, .Machine$double.xmin), 1)        # guard 0/1
      gz  <- sign(gb) * qnorm(gp / 2, lower.tail = FALSE)
      gse <- abs(gb / gz)
      gse[!is.finite(gse) | gz == 0] <- NA_real_
      if (DEBUG && i == 1) message("    (GWAS SE reconstructed from OR + P)")
    } else {
      stop("GWAS has neither SE nor P; cannot obtain standard error")
    }

    # ---- build harmonization frames; effect allele explicit per source ----
    # GWAS: OR/beta signed to A1 (effect = A1).  eQTL: effect = alt_alleles.
    g <- tibble(
      SNP        = as.character(gwas[[g_snp]]),            # keep native case (LD uses lowercase rsIDs)
      CHR        = chr_name,
      BP         = safe_num(gwas[[g_bp]]),
      eff_gwas   = toupper(as.character(gwas[[g_a1]])),   # GWAS effect allele
      noneff_gwas= toupper(as.character(gwas[[g_a2]])),
      beta_gwas  = gb,
      se_gwas    = gse,
      FRQ        = if (!is.na(g_frq)) safe_num(gwas[[g_frq]]) else NA_real_)
    e_effect <- if (!is.na(e_eff)) toupper(as.character(eqtl[[e_eff]]))
                else toupper(as.character(eqtl[[e_a2]]))   # fallback: A2 = effect
    e_other  <- ifelse(e_effect == toupper(as.character(eqtl[[e_a1]])),
                       toupper(as.character(eqtl[[e_a2]])),
                       toupper(as.character(eqtl[[e_a1]])))
    e <- tibble(
      SNP        = as.character(eqtl[[e_snp]]),            # keep native case
      eff_eqtl   = e_effect,                              # eQTL effect allele
      noneff_eqtl= e_other,
      beta_eqtl  = safe_num(eqtl[[e_beta]]),
      se_eqtl    = safe_num(eqtl[[e_se]]),
      FRQ_eqtl   = if (!is.na(e_frq)) safe_num(eqtl[[e_frq]]) else NA_real_)

    overlap <- dplyr::inner_join(g, e, by = "SNP")
    if (debug <- DEBUG) message("    rsID overlap: ", nrow(overlap))
    if (nrow(overlap) < MIN_OVERLAP) { log_target("too_few_overlap_pre_LD", nrow(overlap)); next }

    # ---- pick best LD block by SNP coverage (case-insensitive match) ----
    # LD rownames are lowercase rsIDs; normalize SNP IDs to match.
    overlap$SNP <- tolower(overlap$SNP)
    blk <- ld_blocks %>% filter(chr == chr_name)
    if (nrow(blk) == 0) { log_target("no_LD_block_for_chr", nrow(overlap)); next }
    best <- NULL; best_n <- 0
    for (j in seq_len(nrow(blk))) {
      ld <- tryCatch(readRDS(blk$file[j]), error = function(ex) NULL); if (is.null(ld)) next
      n <- length(intersect(overlap$SNP, tolower(rownames(ld))))
      if (n > best_n) { best_n <- n; best <- j }
    }
    if (is.null(best) || best_n < MIN_OVERLAP) { log_target("LD_block_too_few_snps", best_n); next }

    ld <- readRDS(blk$file[best])
    rownames(ld) <- tolower(rownames(ld)); colnames(ld) <- tolower(colnames(ld))
    shared <- intersect(overlap$SNP, rownames(ld))
    ob <- overlap[overlap$SNP %in% shared, , drop = FALSE]

    # ---- harmonize alleles & betas to LD reference ----
    # MAF for coloc's eQTL dataset: prefer eQTL-sample frequency (FRQ_eqtl),
    # fall back to GWAS FRQ if eQTL lacks one. min(f,1-f) is orientation-invariant.
    frq_src <- ifelse(!is.na(ob$FRQ_eqtl), ob$FRQ_eqtl, ob$FRQ)
    ob$MAF  <- pmin(frq_src, 1 - frq_src)
    ob <- harmonize(ob, ld, debug = DEBUG)
    # drop any SNPs still lacking a usable MAF (coloc rejects NA/0/1 MAF)
    ob <- ob[!is.na(ob$MAF) & ob$MAF > 0 & ob$MAF < 1, , drop = FALSE]
    if (nrow(ob) < MIN_OVERLAP) { log_target("too_few_after_harmonize", nrow(ob)); next }

    # ---- order LD and data identically (single source of truth) ----
    ord <- rownames(ld)[rownames(ld) %in% ob$SNP]
    ob  <- ob[match(ord, ob$SNP), , drop = FALSE]
    ld_sub <- ld[ord, ord, drop = FALSE]

    row <- tryCatch(
      run_one(ob, ld_sub, chr_name, gene_ens, gene_name,
              N_GWAS, N_EQTL, MAX_CS, debug = DEBUG),
      error = function(ex) { message("    run_one error: ", conditionMessage(ex)); NULL })
    if (is.null(row)) { log_target("run_one_returned_null", nrow(ob)); next }

    row$ld_block <- basename(blk$file[best])
    pr <- attr(row, "pairs")
    if (!is.null(pr) && nrow(pr) > 0) { pr$ld_block <- basename(blk$file[best]); all_pairs[[length(all_pairs) + 1]] <- pr }
    all_results[[length(all_results) + 1]] <- row
    if (is_target) log_target(paste0("PROCESSED:", row$status), nrow(ob))
    if (DEBUG) message("    status: ", row$status, " | PP.H4 = ", row$PP_H4)

    # incremental TSV save
    tryCatch(readr::write_tsv(dplyr::bind_rows(all_results), OUT_TSV),
             error = function(ex) message("    TSV write failed: ", conditionMessage(ex)))
  }

  final <- if (length(all_results)) dplyr::bind_rows(all_results) else tibble()
  pairs <- if (length(all_pairs)) dplyr::bind_rows(all_pairs) else
    tibble(gene_name = character(), chr = character(), gwas_cs = integer(),
           eqtl_cs = integer(), gwas_lead = character(), gwas_lead_pip = double(),
           eqtl_lead = character(), eqtl_lead_pip = double(),
           PP_H4 = double())
  tlog  <- if (length(target_log)) dplyr::bind_rows(target_log) else
    tibble(gene = character(), reason = character(), n_snps = integer())

  list(results = final, pairs = pairs, target_log = tlog)
}

# ---------------------------------------------------------------
#  Excel export (multi-sheet)
# ---------------------------------------------------------------
write_excel <- function(out, path) {
  res <- out$results
  wb <- createWorkbook()
  hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                     halign = "center", border = "bottom", fontName = "Arial")
  body <- createStyle(fontName = "Arial")

  add_sheet <- function(name, df) {
    addWorksheet(wb, name)
    if (nrow(df) == 0) { writeData(wb, name, "No rows"); return() }
    writeData(wb, name, df, headerStyle = hdr)
    addStyle(wb, name, body, rows = 2:(nrow(df) + 1),
             cols = 1:ncol(df), gridExpand = TRUE)
    setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
    freezePane(wb, name, firstRow = TRUE)
  }

  add_sheet("All_results", res)
  if (nrow(res) && "PP_H4" %in% names(res))
    add_sheet("Colocalized_H4_gt_0.8", res %>% filter(PP_H4 > 0.8) %>% arrange(desc(PP_H4)))
  # compact summary
  if (nrow(res)) {
    keep <- intersect(c("gene_name","chr","nsnps","status","method","PP_H4","best_H4_snp",
                        "best_H4","gwas_n_cs","eqtl_n_cs","gwas_lambda","eqtl_lambda","ld_block"),
                      names(res))
    add_sheet("Summary", res[, keep, drop = FALSE] %>% arrange(desc(PP_H4)))
  }
  # long-format: one row per colocalizing SuSiE credible-set pair
  # (the standard coloc-SuSiE supplementary-table layout)
  if (!is.null(out$pairs) && nrow(out$pairs) > 0)
    add_sheet("Coloc_pairs", out$pairs %>% arrange(desc(PP_H4)))
  add_sheet("Target_loci_log", out$target_log)

  saveWorkbook(wb, path, overwrite = TRUE)
  message("Wrote ", path)
}

# ---------------------------------------------------------------
#  Diagnostics (runs first; inspect before the full pipeline)
# ---------------------------------------------------------------
diagnose <- function() {
  hr <- function() message(strrep("-", 60))
  message("\n", strrep("=", 60))
  message("  DIAGNOSTIC REPORT")
  message(strrep("=", 60))

  # 1. directories exist
  message("\n[1] Paths")
  message("    INPUT_DIR exists: ", dir.exists(INPUT_DIR), "  (", INPUT_DIR, ")")
  message("    LD_DIR    exists: ", dir.exists(LD_DIR),    "  (", LD_DIR, ")")

  # 2. input files
  gwas_files <- sort(list.files(INPUT_DIR, pattern = "gwas", full.names = TRUE))
  eqtl_files <- sort(list.files(INPUT_DIR, pattern = "eqtl", full.names = TRUE))
  message("\n[2] Input files")
  message("    GWAS files: ", length(gwas_files))
  message("    eQTL files: ", length(eqtl_files))
  if (length(gwas_files) == 0 || length(eqtl_files) == 0) {
    message("    !! No input files found - check INPUT_DIR and 'gwas'/'eqtl' name patterns")
    return(invisible())
  }

  # 3. column names from the first file of each type
  message("\n[3] Column names (first file)")
  g1 <- suppressMessages(readr::read_tsv(gwas_files[1], n_max = 5, show_col_types = FALSE))
  e1 <- suppressMessages(readr::read_tsv(eqtl_files[1], n_max = 5, show_col_types = FALSE))
  message("    GWAS: ", paste(names(g1), collapse = ", "))
  message("    eQTL: ", paste(names(e1), collapse = ", "))
  message("\n    GWAS preview:"); print(utils::head(as.data.frame(g1), 3))
  message("    eQTL preview:"); print(utils::head(as.data.frame(e1), 3))

  # 4. required-column check against the aliases the pipeline uses
  message("\n[4] Required-column resolution")
  req <- list(
    GWAS = list(df = g1, fields = list(
      SNP  = c("SNP","rsid","RSID","variant_id","ID"),
      A1   = c("A1","ALT","EA","effect_allele","Allele1"),
      A2   = c("A2","REF","NEA","other_allele","Allele2"),
      SE   = c("SE","se","StdErr","standard_error","sebeta"),
      BETA = c("BETA","beta","b","Effect","effect","OR","odds_ratio"),
      BP   = c("BP","POS","position","pos","BPos"))),
    eQTL = list(df = e1, fields = list(
      SNP  = c("SNP","rsid","RSID","variant_id","ID"),
      A1   = c("A1","ALT","EA","effect_allele","Allele1"),
      A2   = c("A2","REF","NEA","other_allele","Allele2"),
      SE   = c("SE","se","StdErr","standard_error","sebeta"),
      BETA = c("BETA","beta","b","slope","Effect")))
  )
  for (set in names(req)) {
    message("    ", set, ":")
    for (f in names(req[[set]]$fields)) {
      cands <- req[[set]]$fields[[f]]
      hit   <- cands[cands %in% names(req[[set]]$df)]
      message("      ", formatC(f, width = 5), " -> ",
              if (length(hit)) paste0("FOUND '", hit[1], "'")
              else paste0("MISSING (looked for: ", paste(cands, collapse = "/"), ")"))
    }
  }

  # 5. LD blocks
  message("\n[5] LD blocks")
  ld_blocks <- load_ld_blocks(LD_DIR)
  message("    .rds blocks found: ", nrow(ld_blocks))
  if (nrow(ld_blocks) > 0) {
    message("    chromosomes covered: ",
            paste(sort(unique(ld_blocks$chr)), collapse = ", "))
    ld1 <- tryCatch(readRDS(ld_blocks$file[1]), error = function(e) NULL)
    if (!is.null(ld1)) {
      message("    example block: ", basename(ld_blocks$file[1]),
              "  dim=", nrow(ld1), "x", ncol(ld1))
      message("    example rownames: ", paste(utils::head(rownames(ld1), 3), collapse = ", "), " ...")
    }
  } else {
    message("    !! No LD blocks - every locus will drop with 'no_LD_block_for_chr'")
  }

  # 6. target-locus check
  message("\n[6] Target loci: ", paste(TARGET_LOCI, collapse = ", "))
  for (i in seq_along(gwas_files)) {
    ei <- suppressMessages(readr::read_tsv(eqtl_files[min(i, length(eqtl_files))],
                                           n_max = 5, show_col_types = FALSE))
    gn <- if ("GENESYMBOL" %in% names(ei)) unique(ei$GENESYMBOL)[1] else NA
    if (!is.na(gn) && gn %in% TARGET_LOCI) {
      gi <- suppressMessages(readr::read_tsv(gwas_files[i], show_col_types = FALSE))
      chr_col <- intersect(c("CHR","chr","chromosome"), names(gi))[1]
      chr_name <- if (!is.na(chr_col)) paste0("chr", gsub("^chr","", unique(gi[[chr_col]])[1], ignore.case = TRUE)) else NA
      has_block <- chr_name %in% ld_blocks$chr
      message("    ", gn, " (file ", i, "): chr=", chr_name,
              " | LD block for chr present: ", has_block,
              " | GWAS rows: ", nrow(gi))
    }
  }
  message("\n", strrep("=", 60))
  message("  END DIAGNOSTIC REPORT")
  message("  Set DIAGNOSE_ONLY <- FALSE to run the full pipeline.")
  message(strrep("=", 60), "\n")
  invisible()
}

if (DIAGNOSE_ONLY) {
  diagnose()
  quit(save = "no", status = 0)
}

# ---------------------------------------------------------------
#  Run
# ---------------------------------------------------------------
out <- run_pipeline()
saveRDS(out, OUT_RDS)
write_excel(out, OUT_XLSX)

if (nrow(out$results)) {
  message("\nColoc summary:")
  message("  PP.H4 > 0.8: ", sum(out$results$PP_H4 > 0.8, na.rm = TRUE))
  message("  PP.H4 > 0.5: ", sum(out$results$PP_H4 > 0.5, na.rm = TRUE))
}
if (nrow(out$target_log)) {
  message("\nTarget-locus tracking:")
  print(out$target_log)
}
