# ============================================================
# STEP 2: CLASSIFY SIGNIFICANT GENES INTO
#         Neuronal Communication /
#         Transcriptional Regulation /
#         None/Other
#         based on CC GO annotations
# ============================================================


# ============================================================
# INSTALL / LOAD REQUIRED PACKAGES
# ============================================================

cran_pkgs <- c(
  "dplyr",
  "tidyr",
  "readxl",
  "writexl"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}


# ============================================================
# SETTINGS
# ============================================================

# These are outputs produced by STEP 1
detail_file <- "meta_results_2026_08_25_13_39_56_gene_cluster_assignments.xlsx"
significant_genes_file <- "meta_results_2026_08_25_13_39_56_GO_clusters.xlsx"

# Output from this script
output_file <- "meta_results_2026_08_25_13_39_56_gene_CC_category_classification.xlsx"


# ============================================================
# CHECK INPUT FILE EXISTS
# ============================================================

if (!file.exists(detail_file)) {
  stop(
    "Cannot find input file:\n",
    detail_file,
    "\n\nCurrent working directory is:\n",
    getwd()
  )
}


# ============================================================
# READ GENE-CLUSTER DETAIL SHEET
# ============================================================

gene_detail <- readxl::read_excel(
  detail_file,
  sheet = "Gene_Cluster_Detail"
)

message(
  "Loaded ",
  nrow(gene_detail),
  " gene-GO rows"
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "Gene",
  "GO_Term",
  "Ontology",
  "Cluster"
)

missing_cols <- base::setdiff(
  required_cols,
  colnames(gene_detail)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}


# ============================================================
# GET ALL UNIQUE GENES
# ============================================================

significant_genes <- readxl::read_excel(
  significant_genes_file,
  sheet = "Significant_Genes"
)

if (!"SYMBOL" %in% colnames(significant_genes)) {
  stop("Missing SYMBOL column in Significant_Genes sheet")
}

all_genes <- dplyr::bind_rows(
  significant_genes %>% dplyr::transmute(Gene = as.character(SYMBOL)),
  gene_detail %>% dplyr::transmute(Gene = as.character(Gene))
) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct()

message(
  "Unique genes loaded: ",
  nrow(all_genes)
)


# ============================================================
# KEEP ONLY CELLULAR COMPONENT (CC) ANNOTATIONS
# ============================================================

cc_detail <- gene_detail %>%
  dplyr::filter(
    Ontology == "CC"
  )

message(
  "CC annotation rows: ",
  nrow(cc_detail)
)


# ============================================================
# FLAG RELEVANT GO TERMS PER GENE
#
# GO:0005886 = plasma membrane
# GO:0045202 = synapse
# GO:0005634 = nucleus
# ============================================================

go_flags <- cc_detail %>%
  
  dplyr::group_by(Gene) %>%
  
  dplyr::summarise(
    
    has_plasma_membrane = any(
      GO_Term == "GO:0005886",
      na.rm = TRUE
    ),
    
    has_synapse = any(
      GO_Term == "GO:0045202",
      na.rm = TRUE
    ),
    
    has_nucleus = any(
      GO_Term == "GO:0005634",
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )


# ============================================================
# APPLY CATEGORY RULES
#
# Neuronal Communication:
#   plasma membrane OR synapse
#   AND NOT nucleus
#
# Transcriptional Regulation:
#   nucleus
#   AND NOT plasma membrane
#   AND NOT synapse
#
# Everything else:
#   None/Other
# ============================================================

go_flags <- go_flags %>%
  
  dplyr::mutate(
    
    CC_Category = dplyr::case_when(
      
      (has_plasma_membrane | has_synapse) &
        !has_nucleus ~
        "Neuronal Communication",
      
      has_nucleus &
        !has_plasma_membrane &
        !has_synapse ~
        "Transcriptional Regulation",
      
      TRUE ~
        "None/Other"
    )
  )


# ============================================================
# MERGE CLASSIFICATION BACK TO ALL GENES
#
# Genes without any relevant CC annotation are assigned
# None/Other.
# ============================================================

gene_cc_classified <- all_genes %>%
  
  dplyr::left_join(
    go_flags,
    by = "Gene"
  ) %>%
  
  dplyr::mutate(
    
    has_plasma_membrane = tidyr::replace_na(
      has_plasma_membrane,
      FALSE
    ),
    
    has_synapse = tidyr::replace_na(
      has_synapse,
      FALSE
    ),
    
    has_nucleus = tidyr::replace_na(
      has_nucleus,
      FALSE
    ),
    
    CC_Category = tidyr::replace_na(
      CC_Category,
      "None/Other"
    )
  ) %>%
  
  dplyr::arrange(
    CC_Category,
    Gene
  )


# ============================================================
# CREATE CATEGORY SUMMARY
# ============================================================

cc_category_summary <- gene_cc_classified %>%
  
  dplyr::count(
    CC_Category,
    name = "n_genes"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(n_genes)
  )


# ============================================================
# PRINT SUMMARY
# ============================================================

cat("\n")
cat("============================================\n")
cat("  CC CATEGORY CLASSIFICATION SUMMARY\n")
cat("============================================\n")

print(
  cc_category_summary
)


# ============================================================
# PRINT NEURONAL COMMUNICATION GENES
# ============================================================

cat("\n")
cat("============================================\n")
cat("  NEURONAL COMMUNICATION GENES\n")
cat("============================================\n")

neuronal_genes <- gene_cc_classified %>%
  
  dplyr::filter(
    CC_Category == "Neuronal Communication"
  ) %>%
  
  dplyr::pull(Gene)

print(
  neuronal_genes
)


# ============================================================
# PRINT TRANSCRIPTIONAL REGULATION GENES
# ============================================================

cat("\n")
cat("============================================\n")
cat("  TRANSCRIPTIONAL REGULATION GENES\n")
cat("============================================\n")

transcriptional_genes <- gene_cc_classified %>%
  
  dplyr::filter(
    CC_Category == "Transcriptional Regulation"
  ) %>%
  
  dplyr::pull(Gene)

print(
  transcriptional_genes
)


# ============================================================
# PRINT NONE/OTHER GENES
# ============================================================

cat("\n")
cat("============================================\n")
cat("  NONE/OTHER GENES\n")
cat("============================================\n")

other_genes <- gene_cc_classified %>%
  
  dplyr::filter(
    CC_Category == "None/Other"
  ) %>%
  
  dplyr::pull(Gene)

print(
  other_genes
)


# ============================================================
# SAVE RESULTS
# ============================================================

writexl::write_xlsx(
  
  list(
    
    "CC_Classification" =
      gene_cc_classified,
    
    "CC_Category_Summary" =
      cc_category_summary
    
  ),
  
  output_file
)


# ============================================================
# FINISHED
# ============================================================

message("\n============================================")
message("STEP 2 COMPLETE")
message("============================================")

message(
  "Saved to: ",
  output_file
)

message(
  "Total genes classified: ",
  nrow(gene_cc_classified)
)

message(
  "Neuronal Communication: ",
  sum(
    gene_cc_classified$CC_Category ==
      "Neuronal Communication"
  )
)

message(
  "Transcriptional Regulation: ",
  sum(
    gene_cc_classified$CC_Category ==
      "Transcriptional Regulation"
  )
)

message(
  "None/Other: ",
  sum(
    gene_cc_classified$CC_Category ==
      "None/Other"
  )
)

message("============================================")
