suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(arrow)
  library(locuszoomr)
  library(biomaRt)
  library(EnsDb.Hsapiens.v75)
  library(scales)
})

# ── User inputs ────────────────────────────────────────────────────────────────
gwas_file           <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/MS/munged_noINFOFilter_GRCh37_discovery_metav3.0.meta.txt.gz"
eqtl_file           <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/bmind_fetal/PC10/ExN/cis_nominal/ExN.cis_qtl_pairs.chr6.parquet"
target_gene_symbol  <- "AHI1"
target_ensembl_gene <- "AHI1_ENSG00000135541"
ldlink_token        <- ""
flank_size          <- 2e6
output_pdf          <- "AHI1_locus_plots.pdf"
index_snp           <- "rs13197384"
label_snps          <- c("index", "rs7750586", "rs2746431", "rs11154801",
                         "rs2327586", "rs2614276", "rs4896153", "rs6928977")



gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/PAU/munged_noINFOFilter_GRCh37_PAU_trans_ancestry_Aug2023.txt.gz"
eqtl_file  <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/bmind_dev/PC10/Glia/cis_nominal/Glia.cis_qtl_pairs.chr3.parquet"
target_gene_symbol <- "GPX1"
target_ensembl_gene <- "GPX1_ENSG00000233276"   # used in input eqtl filter
ldlink_token        <- ""
flank_size          <- 2e6
output_pdf          <- "GPX1_locus_plots.pdf"
index_snp           <- "rs940045"
label_snps          <- c("index", "rs6446275", "rs4955430", "rs3448"
                         )



gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/BIP/munged_70INFOFilter_GRCh37_bip2024_multianc_no23andMe.txt.gz"
gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/MDD/munged_70INFOFilter_GRCh37_pgc-mdd2025_no23andMe_div_v3-49-46-01.txt.gz"
eqtl_file  <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/bmind_dev/PC10/NNC/cis_nominal/NNC.cis_qtl_pairs.chr11.parquet"
target_gene_symbol <- "FADS1"
target_ensembl_gene <- "FADS1_ENSG00000149485"  # prefer Ensembl id; can change
ldlink_token        <- ""
flank_size          <- 2e6
output_pdf          <- "FADS1_MDD_locus_plots.pdf"
index_snp           <- "rs174562"
label_snps          <- c("index", "rs174547","rs174548","rs174541","rs174550","rs174560")





gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/Walk/munged_noINFOFilter_GRCh37_WALK_metaGWAS_09032023.txt"
eqtl_file  <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/bmind_fetal/PC10/ExN/cis_nominal/ExN.cis_qtl_pairs.chr16.parquet"
target_gene_symbol <- "RBL2"
target_ensembl_gene <- "RBL2_ENSG00000103479"  # prefer Ensembl id; can change
ldlink_token        <- ""
flank_size          <- 2e6
output_pdf          <- "RBL2_walk_locus_plots.pdf"
index_snp           <- "rs17800727"
label_snps          <- c("index")



# ── 1) Read & filter eQTL ─────────────────────────────────────────────────────
eqtl_raw <- read_parquet(eqtl_file)
eqtl_raw <- eqtl_raw[eqtl_raw$phenotype_id == target_ensembl_gene,
                     c("variant_id", "pval_nominal", "slope", "phenotype_id", "af")]
colnames(eqtl_raw) <- c("SNP", "P", "BETA", "GENE", "AF")

# ── 2) Get SNP positions via GRCh37 BioMart (with https + retry) ──────────────
snp_mart <- useMart(
  biomart = "ENSEMBL_MART_SNP",
  dataset = "hsapiens_snp",
  host    = "https://grch37.ensembl.org"
)

snp_pos <- NULL
for (i in 1:5) {
  tryCatch({
    message("BioMart attempt ", i, " of 5...")
    snp_pos <- getBM(
      attributes = c("refsnp_id", "chr_name", "chrom_start"),
      filters    = "snp_filter",
      values     = eqtl_raw$SNP,
      mart       = snp_mart
    )
    message("BioMart succeeded.")
    break
  }, error = function(e) {
    message("Attempt ", i, " failed: ", e$message)
    if (i < 5) {
      message("Waiting 15 seconds before retry...")
      Sys.sleep(15)
    } else {
      stop("All 5 BioMart attempts failed. Try again later or switch to a local SNP database.")
    }
  })
}

colnames(snp_pos) <- c("SNP", "CHR", "BP")
eqtl <- merge(eqtl_raw, snp_pos, by = "SNP", all.x = TRUE)

# ── 3) Read GWAS ───────────────────────────────────────────────────────────────
gwas  <- read.delim(gwas_file, stringsAsFactors = FALSE)
gwas2 <- gwas[, c("CHR", "BP", "P", "SNP")]

# ── 4) Build locus objects ─────────────────────────────────────────────────────
loc_gwas <- locus(
  data      = gwas2,
  chrom     = "CHR",
  pos       = "BP",
  p         = "P",
  index_snp = index_snp,
  gene      = target_gene_symbol,
  flank     = flank_size,
  ens_db    = EnsDb.Hsapiens.v75
)

eqtl2      <- eqtl[, c("CHR", "BP", "P", "SNP")]
eqtl2$CHR  <- as.character(eqtl2$CHR)
eqtl2$BP   <- as.numeric(eqtl2$BP)

loc_eqtl <- locus(
  data      = eqtl2,
  chrom     = "CHR",
  pos       = "BP",
  p         = "P",
  index_snp = index_snp,
  gene      = target_gene_symbol,
  flank     = flank_size,
  ens_db    = EnsDb.Hsapiens.v75
)

# ── 5) Add LD ──────────────────────────────────────────────────────────────────
message("Fetching LD for GWAS locus...")
loc_gwas_ld <- link_LD(loc_gwas,  token = ldlink_token, method = "proxy")

message("Fetching LD for eQTL locus...")
loc_eqtl_ld <- link_LD(loc_eqtl, token = ldlink_token, method = "proxy")

# ── 6) Shared x-axis range ────────────────────────────────────────────────────
xrange <- range(
  c(loc_gwas_ld$data$BP, loc_eqtl_ld$data$BP),
  na.rm = TRUE
)

# ── 7) Save stacked PDF ────────────────────────────────────────────────────────
pdf(output_pdf, width = 10, height = 14)

split.screen(c(2, 1))

# ── Top: GWAS ──
screen(1)
locus_plot(
  loc_gwas_ld,
  main   = paste("GWAS", target_gene_symbol),
  labels = label_snps,
  xlim   = xrange
)

# ── Bottom: eQTL ──
screen(2)
locus_plot(
  loc_eqtl_ld,
  main   = paste("eQTL", target_gene_symbol),
  labels = label_snps,
  xlim   = xrange
)

close.screen(all = TRUE)
dev.off()

message("Saved: ", output_pdf)

