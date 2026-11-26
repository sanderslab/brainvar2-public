# ============================================================
# 03_motif_summary_plots_HOMER.R
#
# HOMER motif summary plots:
#
# ORIGINAL:
#   1. -log10(p) violin plot
#   2. log2OR violin plot
#   3. significant motif count bar plot
#
# NEW:
#   4. clustered motif x run heatmap - log2OR
#   5. clustered motif x run heatmap - -log10(p)
#   6. clustered motif x run heatmap - binary presence/absence
#   7. run x run Jaccard similarity heatmap
#   8. unique vs shared motif bar plot
#   9. UpSet plot of exact motif-sharing combinations
#
# OUTPUT:
#   PNG + PDF for all plots
#   CSV files for motif sharing and Jaccard statistics
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(tidyr)
library(UpSetR)


# ============================================================
# 1. RESOLVE FUNCTION-NAME CONFLICTS
# ============================================================
#
# Your session contains Bioconductor packages such as IRanges,
# which can conflict with dplyr/base functions.
#
# We explicitly use base::intersect() and base::union() later,
# but also resolve common dplyr conflicts here.
# ============================================================

if ("conflicted" %in% loadedNamespaces()) {
  
  for (fn in c(
    "filter",
    "select",
    "arrange",
    "rename",
    "setdiff",
    "first",
    "last",
    "desc",
    "union",
    "intersect",
    "between",
    "slice"
  )) {
    
    try(
      conflicted::conflicts_prefer(
        getExportedValue("dplyr", fn),
        .quiet = TRUE
      ),
      silent = TRUE
    )
    
  }
  
}


# ============================================================
# 2. CONFIGURATION
# ============================================================

in_file <- "/Users/nrohani/homer_runs/motif_outputs/homer_inputs/HOMER_integrated_results/01_known_motifs_long.csv"

out_dir <- "/Users/nrohani/homer_runs/motif_outputs/homer_inputs/HOMER_integrated_results/motif_plots_out"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# Desired biological ordering
analysis_order <- c(
  "broad_cell_type",
  "development",
  "celltype_development"
)

cell_order <- c(
  "ExN",
  "InN",
  "Glia",
  "NNC"
)

stage_order <- c(
  "prenatal",
  "postnatal"
)



# ============================================================
# MOTIF COLOR PALETTE
# ============================================================

motif_colors <- c(
  "broad_cell_type"      = "#F2A158",  # Cell-type
  "development"          = "#40B272",  # Development
  "celltype_development" = "#40B9ED"   # Cell-type & development
)

motif_labels <- c(
  "broad_cell_type"      = "Cell-type",
  "development"          = "Development",
  "celltype_development" = "Cell-type & development"
)

# Significance/effect-size thresholds
padj_cutoff <- 0.05

log2or_cutoff <- 1

keep_result_type <- "known"


# ============================================================
# 3. READ HOMER DATA
# ============================================================

# Automatically detect delimiter
first_line <- readLines(
  in_file,
  n = 1
)

delim <- if (
  str_detect(first_line, "\t")
) {
  
  "\t"
  
} else if (
  str_detect(first_line, ";")
) {
  
  ";"
  
} else {
  
  ","
  
}


df <- read_delim(
  in_file,
  delim = delim,
  show_col_types = FALSE,
  trim_ws = TRUE
)


# ============================================================
# 4. CHECK REQUIRED COLUMNS
# ============================================================

required_cols <- c(
  "Motif_Name",
  "Consensus",
  "P_value",
  "Log_P_value",
  "q_value_Benjamini",
  "of_Target_Sequences_with_Motif",
  "of_Background_Sequences_with_Motif",
  "result_type",
  "run_name"
)


missing_cols <- base::setdiff(
  required_cols,
  names(df)
)


if (length(missing_cols) > 0) {
  
  stop(
    paste0(
      "Detected delim = '",
      delim,
      "'. Missing expected column(s): ",
      paste(
        missing_cols,
        collapse = ", "
      ),
      "\nColumns actually found: ",
      paste(
        names(df),
        collapse = ", "
      )
    )
  )
  
}


# ============================================================
# 5. HELPER FUNCTION
# ============================================================

# Convert "18.11%" -> 0.1811
pct_to_prop <- function(x) {
  
  as.numeric(
    str_remove(
      x,
      "%"
    )
  ) / 100
  
}


# ============================================================
# 6. CALCULATE OR, log2OR AND -log10(p)
# ============================================================

df <- df %>%
  
  mutate(
    
    # Percentages -> proportions
    pct_target =
      pct_to_prop(
        of_Target_Sequences_with_Motif
      ),
    
    pct_background =
      pct_to_prop(
        of_Background_Sequences_with_Motif
      ),
    
    
    # Prevent Inf when percentage is exactly 0 or 100
    pct_target_adj =
      pmin(
        pmax(
          pct_target,
          1e-4
        ),
        1 - 1e-4
      ),
    
    pct_background_adj =
      pmin(
        pmax(
          pct_background,
          1e-4
        ),
        1 - 1e-4
      ),
    
    
    # Odds ratio
    OR =
      (
        pct_target_adj /
          (1 - pct_target_adj)
      ) /
      (
        pct_background_adj /
          (1 - pct_background_adj)
      ),
    
    
    # log2 odds ratio
    log2OR =
      log2(OR),
    
    
    # HOMER Log_P_value is natural log(p)
    # Convert to -log10(p)
    mlog10P =
      -Log_P_value /
      log(10)
    
  )


# ============================================================
# 7. PARSE run_name
# ============================================================
#
# Handles:
#
# DAR_broad_cell_type_{stage}__{celltype}
#
# DAR_development__{celltype}__{stage}
#
# DAR_celltype_development__{celltype}__{stage}
# ============================================================

df <- df %>%
  
  mutate(
    
    analysis_type = case_when(
      
      str_detect(
        run_name,
        "^DAR_broad_cell_type_"
      ) ~
        "broad_cell_type",
      
      str_detect(
        run_name,
        "^DAR_celltype_development__"
      ) ~
        "celltype_development",
      
      str_detect(
        run_name,
        "^DAR_development__"
      ) ~
        "development",
      
      TRUE ~
        NA_character_
      
    ),
    
    
    stage = if_else(
      
      analysis_type ==
        "broad_cell_type",
      
      str_match(
        run_name,
        "_(prenatal|postnatal)__"
      )[, 2],
      
      str_match(
        run_name,
        "__(prenatal|postnatal)$"
      )[, 2]
      
    ),
    
    
    cell_type = if_else(
      
      analysis_type ==
        "broad_cell_type",
      
      str_match(
        run_name,
        "__(ExN|InN|Glia|NNC)$"
      )[, 2],
      
      str_match(
        run_name,
        "__(ExN|InN|Glia|NNC)__"
      )[, 2]
      
    )
    
  )


if (anyNA(df$analysis_type)) {
  
  warning(
    paste(
      "Some run_name values did not match any expected pattern."
    )
  )
  
}


# ============================================================
# 8. FILTER RESULTS
# ============================================================

if (!is.null(keep_result_type)) {
  
  df <- df[
    df$result_type ==
      keep_result_type,
    ,
    drop = FALSE
  ]
  
}


# Significant + enriched motifs
df <- df[
  !is.na(
    df$q_value_Benjamini
  ) &
    df$q_value_Benjamini <
    padj_cutoff,
  ,
  drop = FALSE
]


df <- df[
  !is.na(
    df$log2OR
  ) &
    df$log2OR >
    log2or_cutoff,
  ,
  drop = FALSE
]


# ============================================================
# 9. DEFINE RUN ORDER
# ============================================================

set_levels <- df %>%
  
  distinct(
    run_name,
    analysis_type,
    stage,
    cell_type
  ) %>%
  
  arrange(
    
    match(
      analysis_type,
      analysis_order
    ),
    
    match(
      stage,
      stage_order
    ),
    
    match(
      cell_type,
      cell_order
    )
    
  ) %>%
  
  pull(
    run_name
  )


# Convert to character for easier clustering
runs <- as.character(
  set_levels
)


# Remove runs not actually present
runs <- runs[
  runs %in%
    unique(
      as.character(
        df$run_name
      )
    )
]


if (length(runs) == 0) {
  
  stop(
    "No runs remain after filtering."
  )
  
}


df <- df %>%
  
  mutate(
    
    run_name =
      factor(
        run_name,
        levels = runs
      ),
    
    analysis_type =
      factor(
        analysis_type,
        levels = analysis_order
      )
    
  )


# ============================================================
# 10. COUNTS OF SIGNIFICANT MOTIFS PER RUN
# ============================================================

n_labels <- df %>%
  
  group_by(
    run_name,
    analysis_type
  ) %>%
  
  summarise(
    
    n_sig = n(),
    
    y_pos_p =
      max(
        mlog10P,
        na.rm = TRUE
      ) * 1.05,
    
    y_pos_lfc =
      max(
        log2OR,
        na.rm = TRUE
      ) * 1.05,
    
    .groups = "drop"
    
  )


# ============================================================
# 11. ORIGINAL PLOT 1:
#     -log10(p) VIOLIN
# ============================================================

p_pval_violin <- ggplot(
  df,
  aes(
    x = run_name,
    y = mlog10P,
    fill = analysis_type
  )
) +
  
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.85,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = motif_colors,
    breaks = analysis_order,
    labels = motif_labels
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 0.6,
    alpha = 0.3,
    color = "black"
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "red",
    position =
      position_dodge(
        width = 0.9
      )
  ) +
  
  geom_text(
    data = n_labels,
    aes(
      x = run_name,
      y = y_pos_p,
      label =
        paste0(
          "n=",
          n_sig
        )
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    vjust = 0
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0.05,
          0.12
        )
      )
  ) +
  
  labs(
    x = "Run",
    y = expression(-log[10](p)),
    title =
      "Motif significance by run",
    subtitle =
      paste0(
        "q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    panel.grid.major.x =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "top"
    
  )


print(
  p_pval_violin
)


# ============================================================
# 12. ORIGINAL PLOT 2:
#     log2OR VIOLIN
# ============================================================

p_lfc <- ggplot(
  df,
  aes(
    x = run_name,
    y = log2OR,
    fill = analysis_type
  )
) +
  
  geom_violin(
    trim = FALSE,
    scale = "width",
    alpha = 0.85,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = motif_colors,
    breaks = analysis_order,
    labels = motif_labels
  ) +
  
  geom_jitter(
    width = 0.15,
    size = 0.6,
    alpha = 0.3,
    color = "black"
  ) +
  
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    linewidth = 0.4,
    color = "red",
    position =
      position_dodge(
        width = 0.9
      )
  ) +
  
  geom_text(
    data = n_labels,
    aes(
      x = run_name,
      y = y_pos_lfc,
      label =
        paste0(
          "n=",
          n_sig
        )
    ),
    inherit.aes = FALSE,
    size = 3,
    fontface = "bold",
    vjust = 0
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0.05,
          0.12
        )
      )
  ) +
  
  labs(
    x = "Run",
    y = expression(log[2](OR)),
    title =
      "Motif fold change by run",
    subtitle =
      paste0(
        "q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    panel.grid.major.x =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "top"
    
  )


print(
  p_lfc
)


# ============================================================
# 13. ORIGINAL PLOT 3:
#     SIGNIFICANT MOTIF COUNTS
# ============================================================

p_bar <- ggplot(
  n_labels,
  aes(
    x = run_name,
    y = n_sig,
    fill = analysis_type
  )
) +
  
  geom_col(
    color = "black",
    linewidth = 0.3,
    alpha = 0.85
  ) +
  scale_fill_manual(
    values = motif_colors,
    breaks = analysis_order,
    labels = motif_labels
  ) +
  
  geom_text(
    aes(
      label = n_sig
    ),
    vjust = -0.4,
    size = 3,
    fontface = "bold"
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0,
          0.12
        )
      )
  ) +
  
  labs(
    x = "Run",
    y = "Number of significant motifs",
    title =
      "Significant motif counts by run",
    subtitle =
      paste0(
        "q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
  ) +
  
  theme_minimal(
    base_size = 12
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    panel.grid.major.x =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      ),
    
    legend.position =
      "top"
    
  )


print(
  p_bar
)


# ============================================================
# 14. SAVE ORIGINAL PLOTS
# ============================================================

ggsave(
  file.path(
    out_dir,
    "mlog10P_violin_by_run.pdf"
  ),
  p_pval_violin,
  width = 10,
  height = 6
)

ggsave(
  file.path(
    out_dir,
    "mlog10P_violin_by_run.png"
  ),
  p_pval_violin,
  width = 10,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    out_dir,
    "log2OR_violin_by_run.pdf"
  ),
  p_lfc,
  width = 10,
  height = 6
)

ggsave(
  file.path(
    out_dir,
    "log2OR_violin_by_run.png"
  ),
  p_lfc,
  width = 10,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    out_dir,
    "n_sig_motifs_bar_by_run.pdf"
  ),
  p_bar,
  width = 10,
  height = 6
)

ggsave(
  file.path(
    out_dir,
    "n_sig_motifs_bar_by_run.png"
  ),
  p_bar,
  width = 10,
  height = 6,
  dpi = 300
)


# ============================================================
# 15. COMBINE ORIGINAL PLOTS
# ============================================================

p_combined <-
  p_pval_violin +
  p_lfc +
  p_bar +
  plot_layout(
    ncol = 3,
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )


print(
  p_combined
)


ggsave(
  file.path(
    out_dir,
    "combined_motif_plots_by_run.pdf"
  ),
  p_combined,
  width = 26,
  height = 10
)

ggsave(
  file.path(
    out_dir,
    "combined_motif_plots_by_run.png"
  ),
  p_combined,
  width = 26,
  height = 6,
  dpi = 300
)


# ============================================================
# ============================================================
# NEW ANALYSIS:
# MOTIF SHARING / OVERLAP
# ============================================================
# ============================================================


# ============================================================
# 16. PREPARE MOTIF x RUN DATA
# ============================================================

motif_run <- df %>%
  
  select(
    Motif_Name,
    Consensus,
    run_name,
    analysis_type,
    stage,
    cell_type,
    log2OR,
    mlog10P
  ) %>%
  
  mutate(
    run_name =
      as.character(
        run_name
      )
  ) %>%
  
  distinct()


# ============================================================
# 17. MOTIF SHARING COUNTS
# ============================================================

motif_sharing <- motif_run %>%
  
  distinct(
    Motif_Name,
    run_name
  ) %>%
  
  count(
    Motif_Name,
    name = "n_runs"
  ) %>%
  
  arrange(
    desc(n_runs),
    Motif_Name
  )


n_total_runs <- length(
  runs
)


motif_sharing <- motif_sharing %>%
  
  mutate(
    
    sharing_class = case_when(
      
      n_runs == 1 ~
        "Unique",
      
      n_runs == n_total_runs ~
        "Shared in all runs",
      
      n_runs > 1 ~
        "Shared",
      
      TRUE ~
        NA_character_
      
    )
    
  )


# ============================================================
# 18. CREATE MOTIF x RUN HEATMAP DATA
# ============================================================

heatmap_df <- motif_run %>%
  
  select(
    Motif_Name,
    run_name,
    log2OR
  ) %>%
  
  distinct() %>%
  
  complete(
    Motif_Name,
    run_name = runs
  )


# ============================================================
# 19. CREATE MATRIX FOR CLUSTERING
# ============================================================

cluster_wide <- heatmap_df %>%
  
  pivot_wider(
    names_from = run_name,
    values_from = log2OR
  )


cluster_df <- as.data.frame(
  cluster_wide
)


rownames(
  cluster_df
) <-
  cluster_df$Motif_Name


cluster_df$Motif_Name <- NULL


# Make sure columns are in biological run order
cluster_df <- cluster_df[
  ,
  runs,
  drop = FALSE
]


# Convert to numeric
cluster_df[] <-
  lapply(
    cluster_df,
    as.numeric
  )


cluster_mat <- as.matrix(
  cluster_df
)


# ============================================================
# 20. CLUSTERING MATRIX
# ============================================================
#
# For clustering:
# NA = motif not significantly enriched
#
# We replace NA by 0 ONLY for calculating distances.
#
# The plotted heatmap retains NA and displays it in grey.
# ============================================================

cluster_mat_for_dist <-
  cluster_mat


cluster_mat_for_dist[
  is.na(
    cluster_mat_for_dist
  )
] <- 0


# Remove completely empty rows
keep_rows <- apply(
  cluster_mat_for_dist,
  1,
  function(x) {
    any(
      is.finite(x)
    )
  }
)


cluster_mat_for_dist <-
  cluster_mat_for_dist[
    keep_rows,
    ,
    drop = FALSE
  ]


# ============================================================
# 21. CLUSTER MOTIF ROWS
# ============================================================

if (
  nrow(
    cluster_mat_for_dist
  ) > 1
) {
  
  row_dist <-
    stats::dist(
      cluster_mat_for_dist,
      method = "euclidean"
    )
  
  row_hclust <-
    stats::hclust(
      row_dist,
      method = "complete"
    )
  
  clustered_motif_order <-
    rownames(
      cluster_mat_for_dist
    )[
      row_hclust$order
    ]
  
} else {
  
  clustered_motif_order <-
    rownames(
      cluster_mat_for_dist
    )
  
}


# ============================================================
# 22. CLUSTER RUN COLUMNS
# ============================================================

if (
  ncol(
    cluster_mat_for_dist
  ) > 1
) {
  
  col_dist <-
    stats::dist(
      t(
        cluster_mat_for_dist
      ),
      method = "euclidean"
    )
  
  col_hclust <-
    stats::hclust(
      col_dist,
      method = "complete"
    )
  
  clustered_run_order <-
    colnames(
      cluster_mat_for_dist
    )[
      col_hclust$order
    ]
  
} else {
  
  clustered_run_order <-
    colnames(
      cluster_mat_for_dist
    )
  
}


# ============================================================
# 23. PRINT CLUSTER ORDER
# ============================================================

cat("\n")
cat("============================================\n")
cat("CLUSTERED RUN ORDER\n")
cat("============================================\n")

print(
  clustered_run_order
)


cat("\n")
cat("============================================\n")
cat("CLUSTERED MOTIF ORDER\n")
cat("============================================\n")

print(
  clustered_motif_order
)


# ============================================================
# 24. APPLY CLUSTERED ORDER TO HEATMAP
# ============================================================

heatmap_df <- heatmap_df %>%
  
  filter(
    
    Motif_Name %in%
      clustered_motif_order,
    
    run_name %in%
      clustered_run_order
    
  ) %>%
  
  mutate(
    
    run_name =
      factor(
        run_name,
        levels =
          clustered_run_order
      ),
    
    Motif_Name =
      factor(
        Motif_Name,
        levels =
          rev(
            clustered_motif_order
          )
      )
    
  )


# ============================================================
# 25. CLUSTERED log2OR HEATMAP
# ============================================================

p_motif_heatmap_log2OR <-
  ggplot(
    heatmap_df,
    aes(
      x = run_name,
      y = Motif_Name,
      fill = log2OR
    )
  ) +
  
  geom_tile(
    color = "white",
    linewidth = 0.1
  ) +
  
  scale_fill_viridis_c(
    option = "magma",
    na.value = "grey90"
  ) +
  
  labs(
    
    x = "Run",
    
    y = "Motif",
    
    fill =
      expression(
        log[2](OR)
      ),
    
    title =
      "Clustered motif enrichment",
    
    subtitle =
      paste0(
        "q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff,
        "; rows and columns hierarchically clustered"
      )
    
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    axis.text.y =
      element_text(
        size = 6
      ),
    
    panel.grid =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
    
  )


print(
  p_motif_heatmap_log2OR
)


# ============================================================
# 26. SAVE log2OR HEATMAP
# ============================================================

ggsave(
  file.path(
    out_dir,
    "motif_enrichment_heatmap_CLUSTERED_log2OR.pdf"
  ),
  p_motif_heatmap_log2OR,
  width = 14,
  height = 12
)


ggsave(
  file.path(
    out_dir,
    "motif_enrichment_heatmap_CLUSTERED_log2OR.png"
  ),
  p_motif_heatmap_log2OR,
  width = 14,
  height = 12,
  dpi = 300
)


# ============================================================
# 27. CREATE -log10(p) HEATMAP DATA
# ============================================================

pval_heatmap_df <- motif_run %>%
  
  select(
    Motif_Name,
    run_name,
    mlog10P
  ) %>%
  
  distinct() %>%
  
  complete(
    Motif_Name,
    run_name = clustered_run_order
  ) %>%
  
  filter(
    Motif_Name %in%
      clustered_motif_order
  ) %>%
  
  mutate(
    
    run_name =
      factor(
        run_name,
        levels =
          clustered_run_order
      ),
    
    Motif_Name =
      factor(
        Motif_Name,
        levels =
          rev(
            clustered_motif_order
          )
      )
    
  )


# ============================================================
# 28. CLUSTERED -log10(p) HEATMAP
# ============================================================

p_motif_heatmap_pval <-
  ggplot(
    pval_heatmap_df,
    aes(
      x = run_name,
      y = Motif_Name,
      fill = mlog10P
    )
  ) +
  
  geom_tile(
    color = "white",
    linewidth = 0.1
  ) +
  
  scale_fill_viridis_c(
    option = "plasma",
    na.value = "grey90"
  ) +
  
  labs(
    
    x = "Run",
    
    y = "Motif",
    
    fill =
      expression(
        -log[10](p)
      ),
    
    title =
      "Clustered motif significance",
    
    subtitle =
      paste0(
        "-log10(p); q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
    
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    axis.text.y =
      element_text(
        size = 6
      ),
    
    panel.grid =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
    
  )


print(
  p_motif_heatmap_pval
)


# ============================================================
# 29. SAVE p-value HEATMAP
# ============================================================

ggsave(
  file.path(
    out_dir,
    "motif_significance_heatmap_CLUSTERED_pvalue.pdf"
  ),
  p_motif_heatmap_pval,
  width = 14,
  height = 12
)


ggsave(
  file.path(
    out_dir,
    "motif_significance_heatmap_CLUSTERED_pvalue.png"
  ),
  p_motif_heatmap_pval,
  width = 14,
  height = 12,
  dpi = 300
)


# ============================================================
# 30. CREATE BINARY PRESENCE/ABSENCE MATRIX
# ============================================================

binary_df <- motif_run %>%
  
  distinct(
    Motif_Name,
    run_name
  ) %>%
  
  mutate(
    present = 1
  ) %>%
  
  complete(
    Motif_Name,
    run_name = clustered_run_order,
    fill = list(
      present = 0
    )
  ) %>%
  
  mutate(
    run_name =
      as.character(
        run_name
      )
  )


# Apply clustered order
binary_df <- binary_df %>%
  
  filter(
    Motif_Name %in%
      clustered_motif_order
  ) %>%
  
  mutate(
    
    run_name =
      factor(
        run_name,
        levels =
          clustered_run_order
      ),
    
    Motif_Name =
      factor(
        Motif_Name,
        levels =
          rev(
            clustered_motif_order
          )
      )
    
  )


# ============================================================
# 31. BINARY PRESENCE/ABSENCE HEATMAP
# ============================================================

p_motif_heatmap_binary <-
  ggplot(
    binary_df,
    aes(
      x = run_name,
      y = Motif_Name,
      fill = factor(
        present
      )
    )
  ) +
  
  geom_tile(
    color = "white",
    linewidth = 0.1
  ) +
  
  scale_fill_manual(
    values = c(
      "0" = "grey90",
      "1" = "black"
    ),
    labels = c(
      "0" = "Absent",
      "1" = "Present"
    ),
    name = "Motif"
  ) +
  
  labs(
    
    x = "Run",
    
    y = "Motif",
    
    title =
      "Clustered motif sharing",
    
    subtitle =
      "Presence/absence of significant enriched motifs"
    
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    axis.text.y =
      element_text(
        size = 6
      ),
    
    panel.grid =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
    
  )


print(
  p_motif_heatmap_binary
)


# ============================================================
# 32. SAVE BINARY HEATMAP
# ============================================================

ggsave(
  file.path(
    out_dir,
    "motif_presence_absence_heatmap_CLUSTERED.pdf"
  ),
  p_motif_heatmap_binary,
  width = 14,
  height = 12
)


ggsave(
  file.path(
    out_dir,
    "motif_presence_absence_heatmap_CLUSTERED.png"
  ),
  p_motif_heatmap_binary,
  width = 14,
  height = 12,
  dpi = 300
)


# ============================================================
# 33. CALCULATE JACCARD SIMILARITY
# ============================================================
#
# IMPORTANT:
# Explicit base::intersect() and base::union() are used.
#
# This prevents the exact conflicted error you encountered.
# ============================================================

motif_sets <- motif_run %>%
  
  distinct(
    run_name,
    Motif_Name
  )


jaccard_df <- expand.grid(
  
  run1 = runs,
  
  run2 = runs,
  
  stringsAsFactors = FALSE
  
) %>%
  
  rowwise() %>%
  
  mutate(
    
    set1 = list(
      motif_sets$Motif_Name[
        motif_sets$run_name ==
          run1
      ]
    ),
    
    set2 = list(
      motif_sets$Motif_Name[
        motif_sets$run_name ==
          run2
      ]
    ),
    
    n_shared =
      length(
        base::intersect(
          set1,
          set2
        )
      ),
    
    n_union =
      length(
        base::union(
          set1,
          set2
        )
      ),
    
    jaccard =
      ifelse(
        n_union == 0,
        NA_real_,
        n_shared /
          n_union
      )
    
  ) %>%
  
  ungroup()


# ============================================================
# 34. JACCARD FACTORS
# ============================================================

jaccard_df <- jaccard_df %>%
  
  mutate(
    
    run1 =
      factor(
        run1,
        levels =
          runs
      ),
    
    run2 =
      factor(
        run2,
        levels =
          runs
      )
    
  )


# ============================================================
# 35. JACCARD HEATMAP
# ============================================================

p_jaccard <-
  ggplot(
    jaccard_df,
    aes(
      x = run1,
      y = run2,
      fill = jaccard
    )
  ) +
  
  geom_tile(
    color = "white",
    linewidth = 0.15
  ) +
  
  geom_text(
    aes(
      label =
        ifelse(
          is.na(jaccard),
          "",
          sprintf(
            "%.2f",
            jaccard
          )
        )
    ),
    size = 2.8
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(
      0,
      1
    ),
    na.value = "grey90"
  ) +
  
  coord_equal() +
  
  labs(
    
    x = "Run",
    
    y = "Run",
    
    fill =
      "Jaccard\nsimilarity",
    
    title =
      "Motif-set similarity between runs",
    
    subtitle =
      paste0(
        "Based on significant enriched motifs; q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
    
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    panel.grid =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
    
  )


print(
  p_jaccard
)


# ============================================================
# 36. SAVE JACCARD HEATMAP
# ============================================================

ggsave(
  file.path(
    out_dir,
    "motif_jaccard_similarity_heatmap.pdf"
  ),
  p_jaccard,
  width = 12,
  height = 10
)


ggsave(
  file.path(
    out_dir,
    "motif_jaccard_similarity_heatmap.png"
  ),
  p_jaccard,
  width = 12,
  height = 10,
  dpi = 300
)


# ============================================================
# 37. UNIQUE VS SHARED MOTIFS
# ============================================================

motif_counts <- motif_run %>%
  
  distinct(
    Motif_Name,
    run_name
  ) %>%
  
  count(
    Motif_Name,
    name = "n_runs"
  )


sharing_by_run <- motif_run %>%
  
  distinct(
    Motif_Name,
    run_name
  ) %>%
  
  left_join(
    motif_counts,
    by = "Motif_Name"
  ) %>%
  
  mutate(
    
    sharing =
      case_when(
        
        n_runs == 1 ~
          "Unique to this run",
        
        n_runs > 1 ~
          "Shared",
        
        TRUE ~
          NA_character_
        
      )
    
  )


sharing_summary <- sharing_by_run %>%
  
  count(
    run_name,
    sharing,
    name = "n_motifs"
  ) %>%
  
  complete(
    
    run_name =
      runs,
    
    sharing = c(
      "Unique to this run",
      "Shared"
    ),
    
    fill = list(
      n_motifs = 0
    )
    
  ) %>%
  
  mutate(
    
    run_name =
      factor(
        run_name,
        levels =
          runs
      )
    
  )


# ============================================================
# 38. UNIQUE VS SHARED BAR PLOT
# ============================================================

p_sharing_bar <-
  ggplot(
    sharing_summary,
    aes(
      x = run_name,
      y = n_motifs,
      fill = sharing
    )
  ) +
  
  geom_col(
    color = "black",
    linewidth = 0.25
  ) +
  
  geom_text(
    aes(
      label =
        ifelse(
          n_motifs == 0,
          "",
          n_motifs
        )
    ),
    position =
      position_stack(
        vjust = 0.5
      ),
    size = 3
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0,
          0.1
        )
      )
  ) +
  
  labs(
    
    x = "Run",
    
    y =
      "Number of significant motifs",
    
    fill =
      "Motif class",
    
    title =
      "Unique and shared motifs by run",
    
    subtitle =
      paste0(
        "q < ",
        padj_cutoff,
        ", log2OR > ",
        log2or_cutoff
      )
    
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
    
    panel.grid.major.x =
      element_blank(),
    
    plot.title =
      element_text(
        face = "bold"
      )
    
  )


print(
  p_sharing_bar
)


# ============================================================
# 39. SAVE UNIQUE/SHARED BAR PLOT
# ============================================================

ggsave(
  file.path(
    out_dir,
    "unique_vs_shared_motifs_by_run.pdf"
  ),
  p_sharing_bar,
  width = 12,
  height = 7
)


ggsave(
  file.path(
    out_dir,
    "unique_vs_shared_motifs_by_run.png"
  ),
  p_sharing_bar,
  width = 12,
  height = 7,
  dpi = 300
)


# ============================================================
# 40. UPSET DATA
# ============================================================

upset_df <- motif_run %>%
  
  distinct(
    Motif_Name,
    run_name
  ) %>%
  
  mutate(
    present = 1
  ) %>%
  
  pivot_wider(
    
    names_from =
      run_name,
    
    values_from =
      present,
    
    values_fill = 0
    
  )


upset_df <-
  as.data.frame(
    upset_df
  )


# Motif names as row names
rownames(
  upset_df
) <-
  upset_df$Motif_Name


# Remove motif column
upset_df$Motif_Name <-
  NULL


# Add any absent run columns
missing_run_columns <-
  base::setdiff(
    runs,
    names(upset_df)
  )


if (
  length(
    missing_run_columns
  ) > 0
) {
  
  for (
    x in
    missing_run_columns
  ) {
    
    upset_df[[x]] <-
      0
    
  }
  
}


# Exact run order
upset_df <-
  upset_df[
    ,
    runs,
    drop = FALSE
  ]


# Make everything numeric
upset_df[] <-
  lapply(
    upset_df,
    as.numeric
  )


# ============================================================
# 41. UPSET PLOT - PDF
# ============================================================

pdf(
  file.path(
    out_dir,
    "motif_sharing_upset_plot.pdf"
  ),
  width = 16,
  height = 10
)


UpSetR::upset(
  
  upset_df,
  
  sets =
    runs,
  
  keep.order =
    TRUE,
  
  nsets =
    length(runs),
  
  nintersects =
    30,
  
  order.by =
    "freq",
  
  decreasing =
    TRUE,
  
  mb.ratio =
    c(
      0.65,
      0.35
    ),
  
  mainbar.y.label =
    "Number of motifs",
  
  sets.x.label =
    "Number of motifs per run"
  
)


dev.off()


# ============================================================
# 42. UPSET PLOT - PNG
# ============================================================

png(
  
  file.path(
    out_dir,
    "motif_sharing_upset_plot.png"
  ),
  
  width = 16,
  height = 10,
  
  units = "in",
  
  res = 300
  
)


UpSetR::upset(
  
  upset_df,
  
  sets =
    runs,
  
  keep.order =
    TRUE,
  
  nsets =
    length(runs),
  
  nintersects =
    30,
  
  order.by =
    "freq",
  
  decreasing =
    TRUE,
  
  mb.ratio =
    c(
      0.65,
      0.35
    ),
  
  mainbar.y.label =
    "Number of motifs",
  
  sets.x.label =
    "Number of motifs per run"
  
)


dev.off()


# ============================================================
# 43. MOTIF SHARING DISTRIBUTION
# ============================================================

sharing_distribution <-
  motif_sharing %>%
  
  count(
    n_runs,
    name = "n_motifs"
  ) %>%
  
  arrange(
    n_runs
  )


# ============================================================
# 44. SAVE SUMMARY TABLES
# ============================================================

write_csv(
  motif_sharing,
  file.path(
    out_dir,
    "motif_sharing_summary.csv"
  )
)


write_csv(
  jaccard_df %>%
    
    select(
      run1,
      run2,
      n_shared,
      n_union,
      jaccard
    ),
  
  file.path(
    out_dir,
    "motif_jaccard_similarity.csv"
  )
)


write_csv(
  sharing_summary,
  file.path(
    out_dir,
    "unique_vs_shared_motif_counts.csv"
  )
)


write_csv(
  sharing_distribution,
  file.path(
    out_dir,
    "motif_sharing_distribution.csv"
  )
)


write_csv(
  
  data.frame(
    cluster_position =
      seq_along(
        clustered_motif_order
      ),
    
    Motif_Name =
      clustered_motif_order
  ),
  
  file.path(
    out_dir,
    "clustered_motif_order.csv"
  )
  
)


write_csv(
  
  data.frame(
    cluster_position =
      seq_along(
        clustered_run_order
      ),
    
    run_name =
      clustered_run_order
  ),
  
  file.path(
    out_dir,
    "clustered_run_order.csv"
  )
  
)


# ============================================================
# 45. FINAL SUMMARY
# ============================================================

cat("\n")
cat("============================================================\n")
cat("MOTIF ANALYSIS COMPLETE\n")
cat("============================================================\n")

cat(
  "Number of runs: ",
  n_total_runs,
  "\n",
  sep = ""
)

cat(
  "Number of motifs: ",
  n_distinct(
    motif_run$Motif_Name
  ),
  "\n",
  sep = ""
)

cat(
  "Unique motifs: ",
  sum(
    motif_sharing$n_runs == 1
  ),
  "\n",
  sep = ""
)

cat(
  "Shared in >=2 runs: ",
  sum(
    motif_sharing$n_runs >= 2
  ),
  "\n",
  sep = ""
)

cat(
  "Shared in ALL runs: ",
  sum(
    motif_sharing$n_runs ==
      n_total_runs
  ),
  "\n",
  sep = ""
)

cat("\n")
cat("Motif-sharing distribution:\n")

print(
  sharing_distribution
)

cat("\n")
cat("Clustered run order:\n")

print(
  clustered_run_order
)

cat("\n")
cat("Output directory:\n")

cat(
  out_dir,
  "\n"
)

cat("============================================================\n")