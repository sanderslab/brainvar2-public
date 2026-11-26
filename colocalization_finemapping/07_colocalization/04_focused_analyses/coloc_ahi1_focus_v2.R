#!/usr/bin/env Rscript
# =====================================================================
#  AHI1 colocalization across ALL cell types x stages -> single Excel.
#  For each of the 8 combinations it:
#    - reconstructs GWAS beta/SE from OR+P
#    - aligns eQTL effect allele (alt_alleles) and detects LD orientation
#      empirically (lower estimate_s_rss lambda wins)
#    - runs coloc.abf (LD-free) and coloc.susie (LD-based)
#    - records both PP.H4 values, the shared lead SNP, CS sizes, lambda
#  Output: AHI1_coloc_all_celltypes.xlsx
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(coloc); library(susieR); library(openxlsx)
})

## ============================ CONFIG ============================
GENE       <- "AHI1"
TRAIT      <- "MS"
CELLTYPES  <- c("ExN", "InN", "Glia", "NNC")
STAGES     <- c("prenatal", "postnatal")
BASE_INPUT <- "/mnt/s3/FXN/QTL/eQTL/coloc_susie"
LD_DIR     <- "/mnt/s3/FXN/wgs_bv1_and_bv2/upload/LD_exact_blocks"
N_GWAS     <- 41505
N_EQTL_BY_STAGE <- c(prenatal = 127, postnatal = 100)
L_SUSIE    <- 5
COVERAGE   <- 0.95
MIN_ABS_CORR <- 0.5
OUT_XLSX   <- "AHI1_coloc_all_celltypes.xlsx"
## ================================================================

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# returns a one-row data.frame for one celltype x stage, or NULL if AHI1 absent
run_combo <- function(CELLTYPE, STAGE) {
  tag <- paste0(CELLTYPE, "_", STAGE)
  cat("\n==================", tag, "==================\n")
  INPUT_DIR <- file.path(BASE_INPUT, paste0(TRAIT, "_bmind_", STAGE, "_", CELLTYPE))
  N_EQTL <- N_EQTL_BY_STAGE[[STAGE]]

  base <- data.frame(
    gene = GENE, cell_type = CELLTYPE, stage = STAGE, n_eqtl = N_EQTL,
    status = NA_character_, n_overlap = NA_integer_, ld_block = NA_character_,
    abf_PP_H3 = NA_real_, abf_PP_H4 = NA_real_,
    susie_PP_H3 = NA_real_, susie_PP_H4 = NA_real_, susie_H4_ratio = NA_real_,
    shared_lead_snp = NA_character_,
    gwas_lead_pip = NA_real_, eqtl_lead_pip = NA_real_,
    gwas_cs_size = NA_integer_, eqtl_cs_size = NA_integer_,
    gwas_lambda = NA_real_, eqtl_lambda = NA_real_,
    ld_orientation = NA_character_, coloc_quality = NA_character_,
    stringsAsFactors = FALSE)

  if (!dir.exists(INPUT_DIR)) { base$status <- "input_dir_missing"; return(base) }

  ef <- sort(list.files(INPUT_DIR, pattern = "eqtl", full.names = TRUE))
  gf <- sort(list.files(INPUT_DIR, pattern = "gwas", full.names = TRUE))
  if (length(ef) == 0 || length(gf) == 0) { base$status <- "no_input_files"; return(base) }

  idx <- which(sapply(ef, function(f) {
    d <- suppressMessages(readr::read_tsv(f, n_max = 5, show_col_types = FALSE))
    "GENESYMBOL" %in% names(d) && GENE %in% d$GENESYMBOL }))[1]
  if (is.na(idx)) { base$status <- "gene_not_found"; return(base) }

  gwas <- suppressMessages(readr::read_tsv(gf[idx], show_col_types = FALSE)) |> as.data.frame()
  eqtl <- suppressMessages(readr::read_tsv(ef[idx], show_col_types = FALSE)) |> as.data.frame()

  # GWAS beta/SE from OR + P
  gb  <- log(safe_num(gwas$OR))
  gp  <- pmin(pmax(safe_num(gwas$P), .Machine$double.xmin), 1)
  gz  <- sign(gb) * qnorm(gp / 2, lower.tail = FALSE)
  gse <- abs(gb / gz); gse[!is.finite(gse)] <- NA_real_
  g <- tibble(SNP = tolower(gwas$SNP), CHR = gwas$CHR, BP = safe_num(gwas$BP),
              beta = gb, se = gse)

  # eQTL effect = alt_alleles
  e <- tibble(SNP = tolower(eqtl$SNP),
              beta = safe_num(eqtl$BETA), se = safe_num(eqtl$SE),
              FRQ = safe_num(eqtl$FRQ))

  m <- inner_join(g, e, by = "SNP", suffix = c("_g", "_e"))
  m <- m[is.finite(m$beta_g) & is.finite(m$se_g) & m$se_g > 0 &
         is.finite(m$beta_e) & is.finite(m$se_e) & m$se_e > 0 &
         is.finite(m$FRQ) & m$FRQ > 0 & m$FRQ < 1, ]
  base$n_overlap <- nrow(m)
  if (nrow(m) < 5) { base$status <- "too_few_overlap"; return(base) }

  # pick LD block
  chr_name <- paste0("chr", gsub("^chr", "", unique(m$CHR)[1], ignore.case = TRUE))
  allb <- list.files(LD_DIR, pattern = "\\.rds$", full.names = TRUE)
  cand <- allb[grep(paste0(chr_name, "[:_]"), basename(allb))]
  if (length(cand) == 0) { base$status <- "no_LD_block"; return(base) }
  best_f <- NULL; best_n <- 0
  for (f in cand) {
    rn <- tolower(rownames(readRDS(f)))
    n <- length(intersect(m$SNP, rn))
    if (n > best_n) { best_n <- n; best_f <- f }
  }
  if (is.null(best_f) || best_n < 5) { base$status <- "LD_too_few_snps"; return(base) }
  base$ld_block <- basename(best_f)

  ld <- readRDS(best_f); rownames(ld) <- tolower(rownames(ld)); colnames(ld) <- rownames(ld)
  shared <- intersect(rownames(ld), m$SNP)
  m  <- m[match(shared, m$SNP), ]
  ld <- ld[shared, shared, drop = FALSE]; ld <- (ld + t(ld)) / 2; diag(ld) <- 1
  ev <- eigen(ld, symmetric = TRUE)
  if (min(ev$values) < 1e-4) {
    ev$values[ev$values < 1e-4] <- 1e-4
    ld <- ev$vectors %*% diag(ev$values) %*% t(ev$vectors)
    d <- sqrt(diag(ld)); ld <- ld / outer(d, d); diag(ld) <- 1
    dimnames(ld) <- list(shared, shared)
  }
  base$n_overlap <- nrow(m)

  # empirical orientation: pick lower-lambda sign per dataset
  pick <- function(beta, se, n) {
    z <- beta / se
    lp <- tryCatch(estimate_s_rss(z =  z, R = ld, n = n), error = function(e) NA_real_)
    ln <- tryCatch(estimate_s_rss(z = -z, R = ld, n = n), error = function(e) NA_real_)
    if (is.na(ln) || (!is.na(lp) && lp <= ln)) list(beta = beta, lambda = lp, flip = FALSE)
    else list(beta = -beta, lambda = ln, flip = TRUE)
  }
  og <- pick(m$beta_g, m$se_g, N_GWAS)
  oe <- pick(m$beta_e, m$se_e, N_EQTL)
  base$gwas_lambda <- round(og$lambda, 4)
  base$eqtl_lambda <- round(oe$lambda, 4)
  base$ld_orientation <- paste0("gwas:", ifelse(og$flip, "flip", "asis"),
                                "|eqtl:", ifelse(oe$flip, "flip", "asis"))
  maf <- pmin(m$FRQ, 1 - m$FRQ)

  # coloc.abf (LD-free)
  abf <- tryCatch(coloc.abf(
    list(beta = og$beta, varbeta = m$se_g^2, snp = m$SNP, type = "quant", N = N_GWAS, sdY = 1),
    list(beta = oe$beta, varbeta = m$se_e^2, snp = m$SNP, type = "quant", N = N_EQTL, MAF = maf)),
    error = function(e) NULL)
  if (!is.null(abf)) {
    base$abf_PP_H3 <- round(as.numeric(abf$summary[["PP.H3.abf"]]), 4)
    base$abf_PP_H4 <- round(as.numeric(abf$summary[["PP.H4.abf"]]), 4)
    if (!is.null(abf$results) && "SNP.PP.H4" %in% names(abf$results)) {
      bi <- which.max(abf$results$SNP.PP.H4)
      base$shared_lead_snp <- as.character(abf$results$snp[bi])
    }
  }

  # coloc.susie (LD-based)
  Dg <- list(beta = og$beta, varbeta = m$se_g^2, snp = m$SNP, position = m$BP,
             type = "quant", N = N_GWAS, sdY = 1, LD = ld)
  De <- list(beta = oe$beta, varbeta = m$se_e^2, snp = m$SNP, position = m$BP,
             type = "quant", N = N_EQTL, MAF = maf, LD = ld)
  sg <- tryCatch(runsusie(Dg, L = L_SUSIE, coverage = COVERAGE, min_abs_corr = MIN_ABS_CORR,
                 check_prior = FALSE, estimate_residual_variance = FALSE,
                 maxit = 1000, repeat_until_convergence = TRUE), error = function(e) NULL)
  se_ <- tryCatch(runsusie(De, L = L_SUSIE, coverage = COVERAGE, min_abs_corr = MIN_ABS_CORR,
                 check_prior = FALSE, estimate_residual_variance = FALSE,
                 maxit = 1000, repeat_until_convergence = TRUE), error = function(e) NULL)
  if (!is.null(sg) && !is.null(se_)) {
    cs <- tryCatch(coloc.susie(sg, se_), error = function(e) NULL)
    if (!is.null(cs) && !is.null(cs$summary) && nrow(cs$summary) > 0) {
      sm <- as.data.frame(cs$summary)
      bp <- which.max(sm$PP.H4.abf)
      base$susie_PP_H3 <- round(sm$PP.H3.abf[bp], 4)
      base$susie_PP_H4 <- round(sm$PP.H4.abf[bp], 4)
      base$susie_H4_ratio <- round(sm$PP.H4.abf[bp] / (sm$PP.H3.abf[bp] + sm$PP.H4.abf[bp]), 4)
      # lead SNPs + CS sizes for the best pair
      lead <- function(s, lab) {
        nm <- names(s$sets$cs); pos <- which(nm %in% c(as.character(lab), paste0("L", lab)))
        if (length(pos) == 0) return(c(snp = NA, pip = NA, size = NA))
        ix <- s$sets$cs[[pos[1]]]; pip <- s$pip[ix]; k <- which.max(pip)
        c(snp = names(s$pip)[ix][k], pip = round(unname(pip[k]), 4), size = length(ix))
      }
      lg <- lead(sg, sm$idx1[bp]); le <- lead(se_, sm$idx2[bp])
      base$gwas_lead_pip <- as.numeric(lg["pip"]); base$eqtl_lead_pip <- as.numeric(le["pip"])
      base$gwas_cs_size  <- as.integer(lg["size"]); base$eqtl_cs_size <- as.integer(le["size"])
      if (is.na(base$shared_lead_snp)) base$shared_lead_snp <- as.character(le["snp"])
    }
  }

  # overall quality call (prefer agreement of the two methods)
  h4 <- max(base$abf_PP_H4 %||% NA, base$susie_PP_H4 %||% NA, na.rm = TRUE)
  lam <- max(base$gwas_lambda %||% NA, base$eqtl_lambda %||% NA, na.rm = TRUE)
  both_high <- (!is.na(base$abf_PP_H4) && base$abf_PP_H4 >= 0.8) &&
               (!is.na(base$susie_PP_H4) && base$susie_PP_H4 >= 0.8)
  base$coloc_quality <- if (is.na(h4)) "none"
    else if (both_high) "high_both_methods"
    else if (!is.na(base$abf_PP_H4) && base$abf_PP_H4 >= 0.8) "high_abf_only"
    else if (h4 >= 0.5) "moderate"
    else "no_coloc"
  base$status <- "ok"
  cat("  abf H4 =", base$abf_PP_H4, "| susie H4 =", base$susie_PP_H4,
      "| lambda =", base$eqtl_lambda, "| quality =", base$coloc_quality, "\n")
  base
}

# ---- run all 8 ----
rows <- list()
for (CT in CELLTYPES) for (ST in STAGES)
  rows[[length(rows) + 1]] <- tryCatch(run_combo(CT, ST),
    error = function(e) { message("combo error ", CT, " ", ST, ": ", e$message); NULL })
res <- dplyr::bind_rows(rows)

# ---- write Excel ----
wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                   halign = "center", border = "bottom", fontName = "Arial")
body <- createStyle(fontName = "Arial")
hi  <- createStyle(fgFill = "#C6EFCE", fontName = "Arial")  # green = high coloc

add_sheet <- function(name, df) {
  addWorksheet(wb, name)
  writeData(wb, name, df, headerStyle = hdr)
  if (nrow(df) > 0) {
    addStyle(wb, name, body, rows = 2:(nrow(df)+1), cols = 1:ncol(df), gridExpand = TRUE)
    setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
    freezePane(wb, name, firstRow = TRUE)
    # highlight high-confidence rows
    if ("coloc_quality" %in% names(df)) {
      hir <- which(df$coloc_quality %in% c("high_both_methods"))
      if (length(hir)) addStyle(wb, name, hi, rows = hir + 1, cols = 1:ncol(df), gridExpand = TRUE)
    }
  }
}

add_sheet("AHI1_all_combos", res)
ok <- res %>% filter(status == "ok") %>%
  arrange(stage, factor(cell_type, levels = CELLTYPES))
add_sheet("AHI1_results_clean", ok)

saveWorkbook(wb, OUT_XLSX, overwrite = TRUE)
cat("\nWrote", OUT_XLSX, "\n")
print(res[, c("cell_type","stage","abf_PP_H4","susie_PP_H4","eqtl_lambda","coloc_quality")])
