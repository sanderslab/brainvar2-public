## --------------------------------------------------
## AHI1 motif disruption
## Variant: chr6:135497759 C>A
## --------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
pkgs <- c(
  "motifbreakR",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BiocParallel",
  "GenomicRanges"
)
for (p in pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    BiocManager::install(p, ask = FALSE, update = FALSE)
  }
}
library(motifbreakR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BiocParallel)
library(GenomicRanges)

chr <- "chr6"
pos <- 135497759L
ref <- "C"
alt <- "A"

## 1) Check the reference base in hg38
base_at_pos <- as.character(
  getSeq(
    BSgenome.Hsapiens.UCSC.hg38,
    GRanges(seqnames = chr, ranges = IRanges(start = pos, end = pos))
  )
)
cat("Reference base at", chr, pos, "in hg38:", base_at_pos, "\n")
if (toupper(base_at_pos) != ref) {
  stop(
    paste0(
      "Reference allele mismatch: genome has ", base_at_pos,
      " but you provided ", ref, ". ",
      "Check genome build, strand, or allele orientation."
    )
  )
}

## 2) Write a clean BED file
## NOTE: motifbreakR requires the "name" field (V4) to be
## chromosome:start:REF:ALT using the 1-based genomic position
## (NOT the 0-based BED start in V2). Strand (V6) does not need to
## match the gene's actual strand -- motifbreakR scores both
## strands of every PWM regardless of what's in this column.
bed_df <- data.frame(
  V1 = chr,
  V2 = pos - 1L,  # BED start is 0-based
  V3 = pos,       # BED end is 1-based
  V4 = paste0(chr, ":", pos, ":", ref, ":", alt),
  V5 = 0,
  V6 = "+"
)
bed_file <- tempfile(fileext = ".bed")
write.table(
  bed_df,
  file = bed_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
cat("BED file written to:", bed_file, "\n")
print(read.delim(bed_file, header = FALSE))

## 3) Import the SNP
snps <- snps.from.file(
  file = bed_file,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  format = "bed"
)
cat("Imported SNP object:\n")
print(snps)
if (length(snps) == 0) {
  stop("snps.from.file returned 0 variants. The BED formatting or allele orientation is wrong.")
}

## 4) Load motifs, restricted to a target TF set
## Start with hocomoco only, broaden to motifbreakR_motif (HOCOMOCO +
## Homer + Factorbook + ENCODE, ~2800 motifs) if hocomoco alone comes
## up empty at a reasonable threshold.
data(hocomoco, package = "motifbreakR")

## TFs directly implicated for AHI1 + the prenatal excitatory-neuron
## GRN TF set. Union'd and de-duped (RFX2 appears in both). HOCOMOCO
## gene symbols are upper-case human symbols, so we upper-case
## everything here to match regardless of how each list was cased.
target_tfs <- toupper(unique(c(
  # AHI1 target-gene TFs
  "CREM", "RFX2", "KLF12",
  # prenatal excitatory-neuron GRN
  "ZBTB18", "ETV6", "MEIS1", "MEIS2", "JDP2", "ZNF708", "FEZF2",
  "PKNOX1", "RFX1", "ELF2", "FOXP2", "BARX2", "STAT1", "ETV2",
  "ZNF140", "NR6A1", "HEY1", "KLF3", "ETV1", "BACH1", "ZNF341",
  "TFAP2E", "MEF2C", "E2F7", "IRF2", "PAX6", "SP3", "E2F6", "LEF1",
  "STAT3", "CREB1", "NFAT5", "TEAD1", "ZNF652", "MEF2A", "IRF7",
  "BATF", "ZBTB11", "ATF2", "PATZ1", "ZNF610", "FOXJ3", "MEIS3",
  "INSM1", "SP1", "ZNF320", "ZNF85", "ZNF184", "EGR4", "ZNF148",
  "ZSCAN16", "DUX4", "KLF15", "RFX2", "ZIC5", "CTCF", "ZNF675",
  "EBF1"
)))

## hocomoco motif names look like "Hsapiens-HOCOMOCO-CREM_f1" -- match
## on the geneSymbol metadata column rather than parsing names.
hocomoco_symbols <- toupper(mcols(hocomoco)$geneSymbol)
pwm <- hocomoco[hocomoco_symbols %in% target_tfs]

cat("Requested", length(target_tfs), "TFs; matched",
    length(unique(hocomoco_symbols[hocomoco_symbols %in% target_tfs])),
    "of them to", length(pwm), "HOCOMOCO motif model(s).\n")
missing_tfs <- setdiff(target_tfs, hocomoco_symbols)
if (length(missing_tfs) > 0) {
  cat("Not found in hocomoco (try motifbreakR_motif for a broader search):\n")
  print(missing_tfs)
}
if (length(pwm) == 0) {
  stop("None of the requested TFs matched a hocomoco motif -- check spelling or try motifbreakR_motif.")
}

## 5) Pass 1: loose settings, just to confirm the pipeline runs and
## something comes back. Treat these results as provisional -- do
## NOT report hits from this pass as real candidates.
res_loose <- motifbreakR(
  snpList = snps,
  pwmList = pwm,
  threshold = 0.75,
  filterp = FALSE,
  method = "ic",
  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  BPPARAM = BiocParallel::SerialParam()
)
cat("Pass 1 (loose) result class:\n")
print(class(res_loose))
cat("Pass 1 (loose) result:\n")
print(res_loose)
if (is.null(res_loose) || length(res_loose) == 0) {
  stop(
    "No motif hits at threshold = 0.75 even with filterp = FALSE. ",
    "Try motifbreakR_motif instead of hocomoco (broader motif set), ",
    "and re-confirm genome build and allele orientation."
  )
}

## 6) Pass 2: tighter, p-value-filtered settings for actual reporting.
## This is what you should use for any hits you plan to interpret.
res <- motifbreakR(
  snpList = snps,
  pwmList = pwm,
  threshold = 0.85,     # default-style reporting threshold
  filterp = TRUE,        # rough p-value pre-filter
  method = "ic",
  bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25),
  BPPARAM = BiocParallel::SerialParam()
)
cat("Pass 2 (filtered) result:\n")
print(res)

if (!is.null(res) && length(res) > 0) {
  ## Calculate actual p-values for the filtered hits (2-step process
  ## required by motifbreakR: filterp=TRUE gives a rough score-based
  ## filter; calculatePvalue() fills in the real Refpvalue/Altpvalue)
  ##
  ## This step (TFMsc2pv) can be very slow at full precision. Two ways
  ## to speed it up, both recommended by the package docs:
  ##   1) granularity: rounds the PWM before computing. 1e-4 is fast
  ##      and fine for triaging; use 1e-6/1e-7 only for a final,
  ##      publication-grade number on a small number of hits.
  ##   2) parallelize: swap SerialParam() for MulticoreParam() -- the
  ##      docs explicitly recommend running this step in parallel.
  n_cores <- max(1L, parallel::detectCores() - 1L)
  res <- calculatePvalue(
    res,
    granularity = 1e-4,
    BPPARAM = if (.Platform$OS.type == "windows") {
      BiocParallel::SerialParam()  # MulticoreParam is unix/mac only
    } else {
      BiocParallel::MulticoreParam(workers = n_cores)
    }
  )
  res_df <- as.data.frame(res, row.names = NULL)
  res_df[] <- lapply(res_df, function(col) {
    if (is.list(col)) {
      vapply(col, function(x) paste(x, collapse = ";"), character(1))
    } else {
      col
    }
  })
  cat("Result dimensions:\n")
  print(dim(res_df))
  print(names(res_df))
  write.csv(res_df, "AHI1_chr6_135497759_motifbreakR_results.csv", row.names = FALSE)
  cat("Filtered results written to AHI1_chr6_135497759_motifbreakR_results.csv\n")
} else {
  cat(
    "No hits survived p-value filtering at threshold = 0.85.\n",
    "The loose pass-1 hits above (threshold 0.75, no p-value filter) ",
    "are still available in res_loose if you want to inspect weaker candidates,\n",
    "but treat those as exploratory, not confirmed.\n"
  )
  res_loose_df <- as.data.frame(res_loose, row.names = NULL)
  res_loose_df[] <- lapply(res_loose_df, function(col) {
    if (is.list(col)) {
      vapply(col, function(x) paste(x, collapse = ";"), character(1))
    } else {
      col
    }
  })
  write.csv(res_loose_df, "AHI1_chr6_135497759_motifbreakR_results_LOOSE.csv", row.names = FALSE)
  cat("Loose/exploratory results written to AHI1_chr6_135497759_motifbreakR_results_LOOSE.csv\n")
}
