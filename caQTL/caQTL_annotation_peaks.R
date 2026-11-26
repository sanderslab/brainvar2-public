#!/usr/bin/env Rscript
#
# Annotate caQTL peak coordinates (single "peak_coordinate" column) with:
#   - feature category (Promoter / 5' UTR / 3' UTR / Exon / Intron / Distal),
#     via ChIPseeker::annotatePeak()
#   - nearest gene symbol
#   - signed distance to nearest gene TSS (bp; negative = upstream,
#     positive = downstream, matches ChIPseeker's own sign convention)
#
# Runs automatically across all combinations of:
#   cell types: ExN, Glia, InN, NNC
#   stages:     Postnatal, Prenatal
#
# Install once:
#   install.packages("BiocManager")
#   BiocManager::install(c("ChIPseeker", "GenomicFeatures", "txdbmaker",
#                           "GenomicRanges", "data.table", "org.Hs.eg.db",
#                           "AnnotationDbi"))
#
# Usage:
#   Rscript annotate_caqtl_peaks.R
#
# Input file: tab-separated, must contain a "peak_coordinate" column with
# peak strings like "chr1:920981-921481" (plus whatever other caQTL columns
# you already have - phenotype_id, variant_id, pval_nominal, etc. - all of
# which are preserved and just get three new columns appended).

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

PROMOTER_UP   <- 2000   # bp upstream of TSS
PROMOTER_DOWN <- 2000   # bp downstream of TSS

CELL_TYPES <- c("ExN", "Glia", "InN", "NNC")
STAGES     <- c("Postnatal", "Prenatal")

# Templates - {celltype} and {stage} get substituted in for each run.
CAQTL_DIR <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/caQTL/sig_caQTL"
INPUT_TEMPLATE  <- file.path(CAQTL_DIR, "{stage}_{celltype}_all_chr_1e-5_annotated.tsv")
OUTPUT_TEMPLATE <- file.path(CAQTL_DIR, "{stage}_{celltype}_all_chr_1e-5_peakanno.tsv")

GTF_PATH <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/Enhancer_discovery/gencode.v41.annotation.gtf"

# ---------------------------------------------------------------------------
# 1. Build (or load a cached) TxDb from the GTF
# ---------------------------------------------------------------------------

build_txdb <- function(gtf_path, use_cache = TRUE) {
  cache_path <- paste0(gtf_path, ".txdb.sqlite")
  
  if (use_cache && file.exists(cache_path) &&
      file.info(cache_path)$mtime >= file.info(gtf_path)$mtime) {
    message("Loading cached TxDb from ", cache_path)
    return(loadDb(cache_path))
  }
  
  message("Building TxDb from GTF (this can take a while the first time)...")
  if (requireNamespace("txdbmaker", quietly = TRUE)) {
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
  } else {
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_path, format = "gtf")
  }
  
  if (use_cache) {
    tryCatch(
      saveDb(txdb, cache_path),
      error = function(e) message("Warning: could not cache TxDb: ", conditionMessage(e))
    )
  }
  
  txdb
}

# ---------------------------------------------------------------------------
# 2. Vectorized peak-string parsing -> GRanges of UNIQUE loci
# ---------------------------------------------------------------------------

parse_peaks_to_granges <- function(peak_strings) {
  # peak_strings: character vector like "chr1:920981-921481"
  m <- regmatches(
    peak_strings,
    regexec("^(chr[[:alnum:]_]+):([0-9]+)-([0-9]+)$", peak_strings)
  )
  
  bad <- vapply(m, length, integer(1)) == 0
  if (any(bad)) {
    stop("Could not parse peak string(s), e.g.: ",
         paste(head(peak_strings[bad], 5), collapse = ", "))
  }
  
  parsed <- do.call(rbind, m)
  chrom  <- parsed[, 2]
  start  <- as.integer(parsed[, 3])
  end    <- as.integer(parsed[, 4])
  
  GRanges(
    seqnames = chrom,
    ranges   = IRanges(start = start, end = end),
    PeakID   = peak_strings
  )
}

# ---------------------------------------------------------------------------
# 3. Annotate unique peaks with ChIPseeker
# ---------------------------------------------------------------------------

annotate_unique_peaks <- function(peaks_gr, txdb) {
  peakAnno <- annotatePeak(
    peaks_gr,
    tssRegion = c(-PROMOTER_UP, PROMOTER_DOWN),
    TxDb      = txdb,
    level     = "gene",
    verbose   = FALSE
  )
  
  anno_df <- as.data.frame(peakAnno)
  
  # anno_df$geneId is whatever was in the GTF's gene_id field - for GENCODE
  # that's a versioned Ensembl gene ID, e.g. "ENSG00000223972.5". ChIPseeker
  # only adds a SYMBOL column if annotatePeak() is called with annoDb=...,
  # which we deliberately don't do here (its automatic ID-type detection is
  # unreliable with versioned Ensembl IDs). Instead, map explicitly:
  ens_clean <- sub("\\.[0-9]+$", "", anno_df$geneId)
  gene_symbol <- suppressMessages(
    mapIds(org.Hs.eg.db, keys = ens_clean, keytype = "ENSEMBL",
           column = "SYMBOL", multiVals = "first")
  )
  anno_df$SYMBOL <- unname(gene_symbol)
  
  # Collapse ChIPseeker's detailed annotation strings down to the requested
  # broad category set: Promoter, 5' UTR, 3' UTR, Exon, Intron, Distal
  # Note: ChIPseeker's "Downstream" (near the 3' end of a gene, outside the
  # promoter window) is folded into "Distal" along with "Distal Intergenic".
  simplify_feature <- function(x) {
    x <- as.character(x)
    out <- rep(NA_character_, length(x))
    out[grepl("^Promoter",         x)] <- "Promoter"
    out[grepl("5' UTR",            x)] <- "5' UTR"
    out[grepl("3' UTR",            x)] <- "3' UTR"
    out[grepl("^Exon",             x)] <- "Exon"
    out[grepl("^Intron",           x)] <- "Intron"
    out[grepl("Downstream",        x)] <- "Distal"
    out[grepl("Distal Intergenic", x)] <- "Distal"
    out[is.na(out)]                    <- "Distal"
    out
  }
  
  data.table(
    PeakID              = anno_df$PeakID,
    feature             = simplify_feature(anno_df$annotation),
    nearest_gene_symbol = anno_df$SYMBOL,
    distance_to_gene    = anno_df$distanceToTSS   # neg = upstream, pos = downstream, 0 = overlapping TSS region
  )
}

# ---------------------------------------------------------------------------
# 4. Per-file annotation (adds columns for peak_coordinate)
# ---------------------------------------------------------------------------

annotate_caqtl_file <- function(input_tsv, output_tsv, txdb) {
  df <- fread(input_tsv, sep = "\t")
  stopifnot("peak_coordinate" %in% names(df))
  
  unique_loci <- unique(df$peak_coordinate)
  peaks_gr <- parse_peaks_to_granges(unique_loci)
  
  annotated <- annotate_unique_peaks(peaks_gr, txdb)
  setkey(annotated, PeakID)
  
  idx <- match(df$peak_coordinate, annotated$PeakID)
  df$peak_feature             <- annotated$feature[idx]
  df$peak_nearest_gene_symbol <- annotated$nearest_gene_symbol[idx]
  df$peak_distance_to_gene    <- annotated$distance_to_gene[idx]
  
  fwrite(df, output_tsv, sep = "\t")
  message("Wrote ", output_tsv)
  
  df
}

# ---------------------------------------------------------------------------
# 4b. Peak-type percentage summary for one annotated file
# ---------------------------------------------------------------------------
# NOTE: percentages are computed over UNIQUE peaks (not QTL rows), since the
# same peak can appear on many rows (one per variant tested against it) and
# counting rows would over-weight peaks with more tested variants.

summarize_feature_pct <- function(df, celltype, stage) {
  peak_level <- unique(df[, .(peak_coordinate, peak_feature)])
  
  counts <- peak_level[, .N, by = peak_feature]
  counts[, pct := 100 * N / sum(N)]
  counts[, `:=`(cell_type = celltype, stage = stage, n_unique_peaks = sum(counts$N))]
  setnames(counts, c("peak_feature", "N"), c("feature", "n_peaks"))
  setcolorder(counts, c("stage", "cell_type", "feature", "n_peaks", "n_unique_peaks", "pct"))
  counts[order(-pct)]
}

# ---------------------------------------------------------------------------
# 5. Loop over cell types x stages
# ---------------------------------------------------------------------------

fill_template <- function(template, celltype, stage) {
  path <- gsub("\\{celltype\\}", celltype, template)
  path <- gsub("\\{stage\\}",    stage,    path)
  path
}

run_all <- function(cell_types = CELL_TYPES, stages = STAGES,
                    input_template = INPUT_TEMPLATE,
                    output_template = OUTPUT_TEMPLATE,
                    summary_path = file.path(CAQTL_DIR, "peak_feature_pct_summary.tsv"),
                    gtf_path = GTF_PATH, use_cache = TRUE) {
  
  txdb <- build_txdb(gtf_path, use_cache = use_cache)
  summaries <- list()
  
  for (ct in cell_types) {
    for (st in stages) {
      input_tsv  <- fill_template(input_template, ct, st)
      output_tsv <- fill_template(output_template, ct, st)
      
      if (!file.exists(input_tsv)) {
        message("Skipping ", ct, " / ", st, " - file not found: ", input_tsv)
        next
      }
      
      message("Annotating ", ct, " / ", st, " ...")
      tryCatch(
        {
          df <- annotate_caqtl_file(input_tsv, output_tsv, txdb)
          summaries[[paste(st, ct)]] <- summarize_feature_pct(df, ct, st)
        },
        error = function(e) {
          message("  FAILED for ", ct, " / ", st, ": ", conditionMessage(e))
        }
      )
    }
  }
  
  if (length(summaries) > 0) {
    combined <- rbindlist(summaries)
    fwrite(combined, summary_path, sep = "\t")
    message("Wrote combined peak-feature percentage summary to ", summary_path)
    print(combined)
  }
  
  invisible(summaries)
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

run_all()