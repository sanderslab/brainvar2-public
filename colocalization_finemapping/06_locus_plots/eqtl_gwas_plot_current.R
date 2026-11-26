#!/usr/bin/env Rscript
# GWAS-eQTL Colocalization Visualization Script
# Produces a faceted 2-row x 4-col plot matching Cell-journal style

# -----------------------
# Configuration
# -----------------------

# --- Dataset 1: MS ---
gwas_file   <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/MS/munged_noINFOFilter_GRCh37_discovery_metav3.0.meta.txt.gz"
focal_snp   <- "rs13197384"
target_gene <- "AHI1_ENSG00000135541"
output_file <- "AHI1_MS_coloc_faceted.pdf"
gene_label  <- "AHI1"
chr_eqtl    <- "chr6"
# 
# --- Dataset 2: PAU ---
gwas_file   <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/PAU/munged_noINFOFilter_GRCh37_PAU_trans_ancestry_Aug2023.txt.gz"
focal_snp   <- "rs940045"
target_gene <- "GPX1_ENSG00000233276"
output_file <- "GPX1_PAU_coloc_faceted.pdf"
gene_label  <- "GPX1"
chr_eqtl    <- "chr3"   # adjust to the chromosome your gene is on

# 

# --- Dataset 3: MDD and BIP---

gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/MDD/munged_70INFOFilter_GRCh37_pgc-mdd2025_no23andMe_div_v3-49-46-01.txt.gz"
gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/BIP/munged_70INFOFilter_GRCh37_bip2024_multianc_no23andMe.txt.gz"

focal_snp   <- "rs174562"
target_gene <- "FADS1_ENSG00000149485"  # prefer Ensembl id; can change
output_file <- "FADS1_BIP_coloc_faceted.pdf"
gene_label  <- "FADS1"
chr_eqtl    <- "chr11"   # adjust to the chromosome your gene is on

# --- Dataset 4: Walk ---


gwas_file <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/coloc/data/GWAS_g37/Walk/munged_noINFOFilter_GRCh37_WALK_metaGWAS_09032023.txt"

focal_snp   <- "rs17800727"
target_gene <- "RBL2_ENSG00000103479"  # prefer Ensembl id; can change
output_file <- "RBL2_walk_coloc_faceted.pdf"
gene_label  <- "RBL2"
chr_eqtl    <- "chr16"   # adjust to the chromosome your gene is on


locus_radius  <- 2e6
eqtl_p_thresh <- 1e-5
gwas_p_thresh <- 5e-8

cell_types   <- c("ExN", "InN", "Glia", "NNC")
stages       <- c("fetal", "dev")
stage_labels <- c(fetal = "Prenatal", dev = "Postnatal")

base_dir <- "/Users/nrohani/Documents/Sanders_lab/BrainVar/QTL/bmind"

# -----------------------
# Packages
# -----------------------
required_packages <- c("data.table", "dplyr", "ggplot2", "scales", "arrow", "ggrepel")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(arrow)
  library(ggrepel)
})

# -----------------------
# Helper: detect column name from candidates
# -----------------------
detect_col <- function(df, choices) {
  found <- intersect(choices, colnames(df))
  if (length(found) == 0) return(NA_character_)
  found[1]
}

# -----------------------
# Load & standardise GWAS
# -----------------------
message("Loading GWAS: ", gwas_file)
gwas <- fread(gwas_file, data.table = FALSE)
message("GWAS columns: ", paste(colnames(gwas), collapse = ", "))

chr_col <- detect_col(gwas, c("CHR", "chr", "chrom", "#CHROM", "CHROM"))
pos_col <- detect_col(gwas, c("BP", "POS", "pos", "position", "base_pair_location"))
p_col   <- detect_col(gwas, c("P", "pval", "PVAL", "p_value", "P.value", "p.value"))
snp_col <- detect_col(gwas, c("SNP", "rsid", "RSID", "ID", "MarkerName", "variant_id"))

message(sprintf("Detected columns -> CHR: %s | BP: %s | P: %s | SNP: %s",
                chr_col, pos_col, p_col, snp_col))

# Robust base-R rename (no tidy-eval, handles already-correct names)
safe_rename <- function(df, old, new) {
  if (!is.na(old) && old %in% colnames(df) && old != new) {
    colnames(df)[colnames(df) == old] <- new
  }
  df
}

gwas <- safe_rename(gwas, chr_col, "CHR")
gwas <- safe_rename(gwas, pos_col, "BP")
gwas <- safe_rename(gwas, snp_col, "SNP")

# Coerce types
gwas$CHR <- as.character(gwas$CHR)
gwas$BP  <- as.integer(gwas$BP)

# Identify the p-value column (whatever it ended up being called after safe_rename)
p_col_final <- intersect(colnames(gwas), c("P_raw", "P", "pval", "PVAL", "p_value", "P.value", "p.value"))
if (length(p_col_final) == 0) stop("No p-value column found in GWAS after renaming.")
p_col_final <- p_col_final[1]
message("Using p-value column: ", p_col_final)
gwas[[p_col_final]] <- as.numeric(gwas[[p_col_final]])

# -----------------------
# Locate focal SNP & define locus
# -----------------------
focal_row <- gwas[gwas$SNP == focal_snp, ]
if (nrow(focal_row) == 0) stop("focal_snp '", focal_snp, "' not found in GWAS file.")

chr0 <- focal_row$CHR[1]
pos0 <- focal_row$BP[1]
message(sprintf("Focal SNP %s: chr%s:%d", focal_snp, chr0, pos0))

gwas_loc <- gwas %>%
  dplyr::filter(CHR == chr0,
                BP  >= pos0 - locus_radius,
                BP  <= pos0 + locus_radius)

# Rename p-value column to P.gwas regardless of what it's currently called
colnames(gwas_loc)[colnames(gwas_loc) == p_col_final] <- "P.gwas"

message(sprintf("GWAS locus SNPs: %d", nrow(gwas_loc)))

# -----------------------
# Collect eQTL data across cell types & stages
# -----------------------
all_data <- list()

for (cell in cell_types) {
  for (stage in stages) {
    
    eqtl_file <- paste0(
      base_dir, "_", stage,
      "/PC10/", cell,
      "/cis_nominal/", cell,
      ".cis_qtl_pairs.", chr_eqtl, ".parquet"
    )
    
    if (!file.exists(eqtl_file)) {
      message("Missing: ", eqtl_file)
      next
    }
    
    eqtl_raw <- tryCatch(read_parquet(eqtl_file), error = function(e) {
      message("Error reading: ", eqtl_file, " — ", e$message)
      NULL
    })
    if (is.null(eqtl_raw)) next
    
    needed <- c("variant_id", "pval_nominal", "phenotype_id")
    if (!all(needed %in% colnames(eqtl_raw))) {
      message("Missing columns in: ", eqtl_file)
      next
    }
    
    eqtl <- eqtl_raw
    colnames(eqtl)[colnames(eqtl) == "variant_id"]   <- "SNP"
    colnames(eqtl)[colnames(eqtl) == "pval_nominal"] <- "P.eqtl"
    colnames(eqtl)[colnames(eqtl) == "phenotype_id"] <- "GENE"
    eqtl <- eqtl[grepl(target_gene, eqtl$GENE, fixed = FALSE), ]
    
    if (nrow(eqtl) == 0) {
      message("No rows for gene '", target_gene, "' in: ", eqtl_file)
      next
    }
    
    shared <- intersect(gwas_loc$SNP, eqtl$SNP)
    if (length(shared) == 0) {
      message("No overlapping SNPs for ", cell, "/", stage)
      next
    }
    
    merged <- dplyr::inner_join(
      eqtl     %>% dplyr::filter(SNP %in% shared),
      gwas_loc %>% dplyr::filter(SNP %in% shared),
      by = "SNP"
    ) %>%
      dplyr::mutate(
        neglog10_eqtl = -log10(P.eqtl),
        neglog10_gwas = -log10(P.gwas),
        sig_eqtl  = P.eqtl < eqtl_p_thresh,
        sig_gwas  = P.gwas < gwas_p_thresh,
        is_focal  = SNP == focal_snp,
        cell_type = cell,
        stage     = stage_labels[stage]
      )
    
    all_data[[paste(cell, stage, sep = "_")]] <- merged
    message(sprintf("  %s / %s: %d SNPs merged", cell, stage, nrow(merged)))
  }
}

if (length(all_data) == 0) stop("No data loaded — check file paths and gene name.")

plot_df <- dplyr::bind_rows(all_data)

# -----------------------
# Point categories
# -----------------------
plot_df <- plot_df %>%
  dplyr::mutate(
    category = case_when(
      is_focal            ~ "Focal SNP",
      sig_eqtl & sig_gwas ~ "Sig in both",
      sig_eqtl            ~ "Sig in eQTL only",
      sig_gwas            ~ "Sig in GWAS only",
      TRUE                ~ "Not significant"
    ),
    category  = factor(category, levels = c(
      "Not significant", "Sig in GWAS only", "Sig in eQTL only",
      "Sig in both", "Focal SNP"
    )),
    stage     = factor(stage,     levels = c("Prenatal", "Postnatal")),
    cell_type = factor(cell_type, levels = c("ExN", "InN", "Glia", "NNC"))
  )

cat_colors <- c(
  "Not significant"  = "grey60",
  "Sig in GWAS only" = "#E8392A",
  "Sig in eQTL only" = "#4B9CD3",
  "Sig in both"      = "#F5A623",
  "Focal SNP"        = "#6B4FBB"
)

cat_sizes <- c(
  "Not significant"  = 1.2,
  "Sig in GWAS only" = 1.8,
  "Sig in eQTL only" = 1.8,
  "Sig in both"      = 2.0,
  "Focal SNP"        = 3.5
)

cat_shapes <- c(
  "Not significant"  = 16,
  "Sig in GWAS only" = 16,
  "Sig in eQTL only" = 16,
  "Sig in both"      = 16,
  "Focal SNP"        = 21
)

x_max <- ceiling(max(plot_df$neglog10_eqtl, na.rm = TRUE))
y_max <- ceiling(max(plot_df$neglog10_gwas, na.rm = TRUE))
axis_max <- max(x_max, y_max)   # kept for label positioning only
message(sprintf("Axis limits -> x: 0-%d  y: 0-%d", x_max, y_max))

label_df <- plot_df %>%
  dplyr::filter(is_focal) %>%
  dplyr::select(SNP, neglog10_eqtl, neglog10_gwas, cell_type, stage) %>%
  dplyr::distinct()

# -----------------------
# Plot
# -----------------------
p <- ggplot(plot_df %>% dplyr::arrange(category),
            aes(x = neglog10_eqtl, y = neglog10_gwas)) +
  
  geom_point(aes(color = category, size = category, shape = category),
             alpha = 0.75, stroke = 0) +
  
  geom_point(data = dplyr::filter(plot_df, is_focal),
             shape = 21, size = 3.8, stroke = 1.4,
             color = "#6B4FBB", fill = "#6B4FBB") +
  
  geom_text(data = label_df,
            aes(label = paste0("italic('", SNP, "')")),
            x = 0.5, y = axis_max * 0.97,
            hjust = 0, vjust = 1,
            color = "#6B4FBB", size = 2.8,
            parse = TRUE) +
  
  geom_vline(xintercept = -log10(eqtl_p_thresh),
             linetype = "dashed", linewidth = 0.4, color = "grey40") +
  geom_hline(yintercept = -log10(gwas_p_thresh),
             linetype = "dashed", linewidth = 0.4, color = "grey40") +
  
  facet_grid(stage ~ cell_type) +
  
  scale_color_manual(values = cat_colors, name = NULL,
                     guide = guide_legend(override.aes = list(size = 3, shape = 16))) +
  scale_size_manual( values = cat_sizes,  name = NULL) +
  scale_shape_manual(values = cat_shapes, name = NULL) +
  
  coord_cartesian(xlim = c(0, x_max), ylim = c(0, y_max)) +
  scale_x_continuous(breaks = pretty(c(0, x_max), n = 5)) +
  scale_y_continuous(breaks = pretty(c(0, y_max), n = 5)) +
  
  labs(
    title = bquote(italic(.(gene_label))),
    x     = expression(-log[10](eQTL~italic(P))),
    y     = expression(-log[10](GWAS~italic(P)))
  ) +
  
  theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "italic", size = 12,
                                    margin = margin(b = 6)),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text       = element_text(size = 9, color = "black"),
    strip.text.x     = element_text(margin = margin(b = 3, t = 3)),
    strip.text.y     = element_text(margin = margin(l = 3, r = 3)),
    panel.background = element_rect(fill = "#EBEBEB", color = "grey80"),
    panel.grid.major = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.spacing    = unit(0.35, "lines"),
    axis.text        = element_text(size = 7.5, color = "black"),
    axis.title       = element_text(size = 9),
    axis.ticks       = element_line(linewidth = 0.3, color = "grey50"),
    legend.position  = "right",
    legend.key       = element_rect(fill = NA, color = NA),
    legend.key.size  = unit(0.55, "lines"),
    legend.text      = element_text(size = 8),
    legend.margin    = margin(l = 4),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin      = margin(8, 8, 6, 6)
  ) +
  
  guides(
    color = guide_legend(override.aes = list(size = 3, shape = 16, stroke = 0)),
    size  = "none",
    shape = "none"
  )

# -----------------------
# Save
# -----------------------
fwrite(plot_df, sub("\\.pdf$", "_plot_data.csv.gz", output_file))

ggsave(output_file, p, width = 10, height = 5.5)
message("Saved: ", output_file)