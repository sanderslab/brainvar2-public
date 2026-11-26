# ============================================================
# STEP 1: GO TERM CLUSTERING FOR SIGNIFICANT META-ANALYSIS GENES
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc_pkgs <- c(
  "clusterProfiler",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "GO.db",
  "rrvgo"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

cran_pkgs <- c(
  "dplyr",
  "stringr",
  "readr",
  "writexl"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)
  ) install.packages(pkg)
}

library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(rrvgo)
library(dplyr)
library(stringr)
library(readr)
library(writexl)


# ============================================================
# SETTINGS
# ============================================================

csv_file <- "/Users/nrohani/Downloads/meta_results_2026_08_25_13_39_56.csv"

# Significance cutoff
p_cutoff <- 0.05

# rrvgo similarity threshold
rrvgo_threshold <- 0.7


# ============================================================
# READ META-ANALYSIS RESULTS
# ============================================================

df_input <- read_csv(
  csv_file,
  show_col_types = FALSE
)

message("Loaded ", nrow(df_input), " rows from: ", csv_file)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "Gene",
  "Case-Control SCHEMA P-value"
)

missing_cols <- setdiff(required_cols, colnames(df_input))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ============================================================
# GET SIGNIFICANT GENES ONLY
# ============================================================

df_sig <- df_input %>%
  mutate(
    SCHEMA_P = suppressWarnings(
      as.numeric(`Case-Control SCHEMA P-value`)
    )
  ) %>%
  filter(
    !is.na(SCHEMA_P),
    SCHEMA_P < p_cutoff
  )

message(
  "Significant genes at p < ",
  p_cutoff,
  ": ",
  nrow(df_sig)
)


# ============================================================
# CLEAN ENSEMBL GENE IDS
# ============================================================

gene_list_ensembl <- df_sig$Gene %>%
  as.character() %>%
  str_trim() %>%
  str_replace("\\..*$", "") %>%   # remove version suffix if present
  na.omit() %>%
  unique()

message(
  "Unique significant Ensembl genes: ",
  length(gene_list_ensembl)
)


# ============================================================
# CONVERT ENSEMBL IDS -> HGNC SYMBOL
# ============================================================

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_list_ensembl,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
) %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

message(
  "Significant Ensembl genes mapped to symbols: ",
  nrow(gene_map)
)

unmapped_genes <- setdiff(
  gene_list_ensembl,
  gene_map$ENSEMBL
)

message(
  "Unmapped Ensembl genes: ",
  length(unmapped_genes)
)


# ============================================================
# CREATE SIGNIFICANT GENE TABLE WITH SYMBOLS
# ============================================================

sig_gene_table <- df_sig %>%
  mutate(
    ENSEMBL = Gene %>%
      as.character() %>%
      str_trim() %>%
      str_replace("\\..*$", "")
  ) %>%
  left_join(
    gene_map,
    by = "ENSEMBL"
  )

# Only genes that can be mapped to HGNC symbols
sig_gene_table_mapped <- sig_gene_table %>%
  filter(
    !is.na(SYMBOL),
    SYMBOL != ""
  ) %>%
  distinct(SYMBOL, .keep_all = TRUE)

gene_list <- sig_gene_table_mapped$SYMBOL %>%
  as.character() %>%
  unique()

message(
  "Unique significant genes with valid symbols: ",
  length(gene_list)
)


# ============================================================
# GET ALL GO ANNOTATIONS
# ============================================================

go_all_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = gene_list,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "GO", "ONTOLOGY")
) %>%
  filter(
    !is.na(GO),
    ONTOLOGY %in% c("BP", "CC", "MF")
  )

message(
  "Total GO annotations retrieved: ",
  nrow(go_all_annot)
)

message(
  "  BP: ",
  sum(go_all_annot$ONTOLOGY == "BP"),
  " | CC: ",
  sum(go_all_annot$ONTOLOGY == "CC"),
  " | MF: ",
  sum(go_all_annot$ONTOLOGY == "MF")
)


# ============================================================
# HELPER: RUN rrvgo FOR ONE ONTOLOGY
# ============================================================

run_rrvgo <- function(
    annot_df,
    ontology,
    threshold = 0.7
) {
  
  go_terms <- unique(annot_df$GO)
  
  message(
    "\n[",
    ontology,
    "] ",
    length(go_terms),
    " unique terms - computing similarity matrix..."
  )
  
  # Need at least 2 GO terms
  if (length(go_terms) < 2) {
    message(
      "[",
      ontology,
      "] Fewer than 2 GO terms. Skipping."
    )
    
    return(
      data.frame(
        go = character(),
        term = character(),
        parentTerm = character(),
        Ontology = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  sim_mat <- calculateSimMatrix(
    go_terms,
    orgdb = "org.Hs.eg.db",
    ont = ontology,
    method = "Rel"
  )
  
  scores <- setNames(
    rep(1, length(go_terms)),
    go_terms
  )
  
  reduced <- reduceSimMatrix(
    sim_mat,
    scores,
    threshold = threshold,
    orgdb = "org.Hs.eg.db"
  )
  
  reduced$Ontology <- ontology
  
  reduced
}


# ============================================================
# RUN FOR EACH ONTOLOGY
# ============================================================

bp_annot <- go_all_annot %>%
  filter(ONTOLOGY == "BP")

cc_annot <- go_all_annot %>%
  filter(ONTOLOGY == "CC")

mf_annot <- go_all_annot %>%
  filter(ONTOLOGY == "MF")


reduced_bp <- run_rrvgo(
  bp_annot,
  "BP",
  threshold = rrvgo_threshold
)

reduced_cc <- run_rrvgo(
  cc_annot,
  "CC",
  threshold = rrvgo_threshold
)

reduced_mf <- run_rrvgo(
  mf_annot,
  "MF",
  threshold = rrvgo_threshold
)


# Combine all reduced terms
reduced_all <- bind_rows(
  reduced_bp,
  reduced_cc,
  reduced_mf
)


# ============================================================
# SUMMARISE CLUSTERS PER ONTOLOGY
# ============================================================

summarise_clusters <- function(
    reduced_df,
    ontology
) {
  
  if (nrow(reduced_df) == 0) {
    return(
      data.frame(
        parentTerm = character(),
        Ontology = character(),
        n_terms = integer(),
        example_terms = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  reduced_df %>%
    filter(Ontology == ontology) %>%
    group_by(parentTerm) %>%
    summarise(
      Ontology = ontology,
      n_terms = n(),
      example_terms = paste(
        head(term, 3),
        collapse = " | "
      ),
      .groups = "drop"
    ) %>%
    arrange(desc(n_terms))
}


clusters_bp <- summarise_clusters(
  reduced_all,
  "BP"
)

clusters_cc <- summarise_clusters(
  reduced_all,
  "CC"
)

clusters_mf <- summarise_clusters(
  reduced_all,
  "MF"
)


# ============================================================
# PRINT CLUSTER NAMES
# ============================================================

cat("\n============================================\n")
cat("  BP CLUSTERS (Biological Process)\n")
cat("============================================\n")
print(clusters_bp, n = 60)

cat("\n============================================\n")
cat("  CC CLUSTERS (Cellular Component)\n")
cat("============================================\n")
print(clusters_cc, n = 60)

cat("\n============================================\n")
cat("  MF CLUSTERS (Molecular Function)\n")
cat("============================================\n")
print(clusters_mf, n = 60)


# ============================================================
# SAVE CLUSTER SUMMARY
# ============================================================

write_xlsx(
  list(
    "Significant_Genes" = sig_gene_table_mapped,
    "Unmapped_Genes" = data.frame(
      ENSEMBL = unmapped_genes
    ),
    "BP_clusters" = clusters_bp,
    "CC_clusters" = clusters_cc,
    "MF_clusters" = clusters_mf
  ),
  "meta_results_2026_08_25_13_39_56_GO_clusters.xlsx"
)

message(
  "\nCluster names saved to: ",
  "meta_results_2026_08_25_13_39_56_GO_clusters.xlsx"
)


# ============================================================
# SAVE GENE-LEVEL CLUSTER ASSIGNMENTS
# ============================================================

gene_cluster_detail <- go_all_annot %>%
  
  dplyr::left_join(
    reduced_bp %>%
      dplyr::select(go, parentTerm) %>%
      dplyr::rename(BP_parent = parentTerm),
    by = c("GO" = "go")
  ) %>%
  
  dplyr::left_join(
    reduced_cc %>%
      dplyr::select(go, parentTerm) %>%
      dplyr::rename(CC_parent = parentTerm),
    by = c("GO" = "go")
  ) %>%
  
  dplyr::left_join(
    reduced_mf %>%
      dplyr::select(go, parentTerm) %>%
      dplyr::rename(MF_parent = parentTerm),
    by = c("GO" = "go")
  ) %>%
  
  dplyr::mutate(
    Cluster = dplyr::case_when(
      ONTOLOGY == "BP" & !is.na(BP_parent) ~ BP_parent,
      ONTOLOGY == "CC" & !is.na(CC_parent) ~ CC_parent,
      ONTOLOGY == "MF" & !is.na(MF_parent) ~ MF_parent,
      TRUE ~ NA_character_
    )
  ) %>%
  
  dplyr::select(
    Gene = SYMBOL,
    GO_Term = GO,
    Ontology = ONTOLOGY,
    Cluster
  ) %>%
  
  dplyr::filter(
    !is.na(Cluster)
  ) %>%
  
  dplyr::arrange(
    Gene,
    Ontology,
    Cluster
  )

message(
  "Gene-level cluster assignments created: ",
  nrow(gene_cluster_detail),
  " rows"
)

# ============================================================
# LONG: ONE ROW PER GENE x CLUSTER
# ============================================================

gene_cluster_summary <- gene_cluster_detail %>%
  
  group_by(
    Gene,
    Ontology,
    Cluster
  ) %>%
  
  summarise(
    n_GO_terms = n(),
    .groups = "drop"
  ) %>%
  
  arrange(
    Gene,
    Ontology,
    desc(n_GO_terms)
  )


# ============================================================
# WIDE: ONE ROW PER GENE
# ============================================================

gene_cluster_wide <- gene_cluster_summary %>%
  
  group_by(Gene) %>%
  
  summarise(
    BP_clusters = paste(
      Cluster[Ontology == "BP"],
      collapse = " | "
    ),
    
    CC_clusters = paste(
      Cluster[Ontology == "CC"],
      collapse = " | "
    ),
    
    MF_clusters = paste(
      Cluster[Ontology == "MF"],
      collapse = " | "
    ),
    
    All_clusters = paste(
      unique(Cluster),
      collapse = " | "
    ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    across(
      ends_with("_clusters"),
      ~ ifelse(. == "", NA_character_, .)
    )
  )


# ============================================================
# SAVE GENE-LEVEL OUTPUT
# ============================================================

write_xlsx(
  list(
    "Gene_Cluster_Wide" = gene_cluster_wide,
    "Gene_Cluster_Long" = gene_cluster_summary,
    "Gene_Cluster_Detail" = gene_cluster_detail
  ),
  "meta_results_2026_08_25_13_39_56_gene_cluster_assignments.xlsx"
)

message(
  "\nGene-level cluster assignments saved to: ",
  "meta_results_2026_08_25_13_39_56_gene_cluster_assignments.xlsx"
)

message(
  "  Significant genes: ",
  length(gene_list)
)

message(
  "\nNext step: run the CC category classification script."
)