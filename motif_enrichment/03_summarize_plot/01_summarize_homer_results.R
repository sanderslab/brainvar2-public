# ============================================================
# MASTER HOMER MOTIF ANALYSIS PIPELINE
# ============================================================
#
# INPUT:
#   homer_inputs/
#       *_HOMER/
#           knownResults.txt
#           knownResults.html
#           homerResults.html
#           homerResults/
#           homerMotifs.all.motifs
#           homerMotifs.motifs8
#           homerMotifs.motifs10
#           homerMotifs.motifs12
#           motifFindingParameters.txt
#
# OUTPUT:
#   HOMER_integrated_results/
#
# Main outputs:
#   01_known_motifs_long.csv
#   01_known_motifs_standardized_long.csv
#   02_denovo_motifs_long.csv
#   03_run_metadata.csv
#   04_parameters_long.csv
#   05_motif_recurrence.csv
#   06_motif_stage_shift.csv
#   07_celltype_specificity.csv
#   08_denovo_summary.csv
#   10_summary_by_run.csv
#   11_motif_heatmap_matrix.csv
#   11_motif_enrichment_heatmap.pdf
#   12_motif_rank_stability_long.csv
#   ALL_HOMER_RESULTS.xlsx
#
# ============================================================


# ============================================================
# 0. PACKAGES
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "readr",
  "rvest",
  "xml2",
  "tibble",
  "openxlsx",
  "pheatmap"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(rvest)
  library(xml2)
  library(tibble)
  library(openxlsx)
  library(pheatmap)
  
})


# ============================================================
# 1. USER SETTINGS
# ============================================================

# CHANGE THIS ONLY
homer_root <- "/Users/nrohani/homer_runs/motif_outputs/homer_inputs"


# Output folder
analysis_root <- file.path(
  homer_root,
  "HOMER_integrated_results"
)

dir.create(
  analysis_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# Thresholds
FDR_THRESHOLD <- 0.05

# This is a derived ratio:
# target percentage / background percentage
MIN_TARGET_BACKGROUND_RATIO <- 1.5

# Minimum number of independent HOMER runs
# in which a motif is significant
MIN_RECURRENCE <- 2


# Optional motif-instance analysis
#
# I recommend FALSE initially.
# Turn TRUE later after the core results are validated.
RUN_MOTIF_INSTANCES <- FALSE

HOMER_ANNOTATE <- "annotatePeaks.pl"
HOMER_GENOME <- "hg38"

instance_dir <- file.path(
  analysis_root,
  "motif_instances"
)

dir.create(
  instance_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. BASIC HELPERS
# ============================================================

safe_numeric <- function(x) {
  
  x <- as.character(x)
  
  x <- trimws(x)
  
  x <- gsub(
    ",",
    "",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "%",
    "",
    x,
    fixed = TRUE
  )
  
  x <- gsub(
    "^<",
    "",
    x
  )
  
  x <- gsub(
    "^>",
    "",
    x
  )
  
  suppressWarnings(
    as.numeric(x)
  )
}


clean_colnames_unique <- function(x) {
  
  x <- as.character(x)
  
  x <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    x
  )
  
  x <- gsub(
    "_+",
    "_",
    x
  )
  
  x <- gsub(
    "^_|_$",
    "",
    x
  )
  
  make.unique(
    x,
    sep = "_"
  )
}


# Remove only columns that are 100% NA
drop_all_na_columns <- function(df) {
  
  if (
    is.null(df) ||
    ncol(df) == 0
  ) {
    return(df)
  }
  
  keep <- vapply(
    df,
    function(x) {
      
      any(
        !is.na(x)
      )
      
    },
    logical(1)
  )
  
  df[
    ,
    keep,
    drop = FALSE
  ]
}


# Safely get minimum finite value
safe_min <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  min(
    x,
    na.rm = TRUE
  )
}


# Safely get maximum finite value
safe_max <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  max(
    x,
    na.rm = TRUE
  )
}


# Safely median
safe_median <- function(x) {
  
  x <- x[
    is.finite(x)
  ]
  
  if (length(x) == 0) {
    return(NA_real_)
  }
  
  median(
    x,
    na.rm = TRUE
  )
}


# ============================================================
# 3. PARSE HOMER FOLDER NAME
# ============================================================

parse_homer_metadata <- function(folder_name) {
  
  x <- stringr::str_remove(
    folder_name,
    "_HOMER$"
  )
  
  parts <- stringr::str_split(
    x,
    stringr::fixed("__")
  )[[1]]
  
  tibble::tibble(
    
    dataset =
      if (length(parts) >= 1)
        parts[1]
    else
      NA_character_,
    
    analysis =
      if (length(parts) >= 2)
        parts[2]
    else
      NA_character_,
    
    cell_type =
      if (length(parts) >= 3)
        parts[3]
    else
      NA_character_,
    
    stage =
      if (length(parts) >= 4)
        parts[4]
    else
      NA_character_,
    
    run_name =
      folder_name
  )
}


# ============================================================
# 4. FIND HOMER FOLDERS
# ============================================================

if (!dir.exists(homer_root)) {
  
  stop(
    "\nHOMER directory does not exist:\n",
    homer_root,
    "\n"
  )
}


homer_dirs <- list.dirs(
  homer_root,
  recursive = FALSE,
  full.names = TRUE
)

homer_dirs <- homer_dirs[
  stringr::str_detect(
    basename(homer_dirs),
    "_HOMER$"
  )
]


if (length(homer_dirs) == 0) {
  
  stop(
    "\nNo *_HOMER directories were found in:\n",
    homer_root,
    "\n"
  )
}


cat(
  "\n============================================\n"
)

cat(
  "Found ",
  length(homer_dirs),
  " HOMER result folders\n",
  sep = ""
)

cat(
  "============================================\n\n"
)


# ============================================================
# 5. FILE INVENTORY
# ============================================================

run_metadata <- purrr::map_dfr(
  homer_dirs,
  function(folder) {
    
    meta <- parse_homer_metadata(
      basename(folder)
    )
    
    tibble::tibble(
      
      source_folder =
        folder,
      
      knownResults_txt =
        file.exists(
          file.path(
            folder,
            "knownResults.txt"
          )
        ),
      
      knownResults_html =
        file.exists(
          file.path(
            folder,
            "knownResults.html"
          )
        ),
      
      homerResults_html =
        file.exists(
          file.path(
            folder,
            "homerResults.html"
          )
        ),
      
      homerResults_directory =
        dir.exists(
          file.path(
            folder,
            "homerResults"
          )
        ),
      
      homerMotifs_all =
        file.exists(
          file.path(
            folder,
            "homerMotifs.all.motifs"
          )
        ),
      
      homerMotifs_8 =
        file.exists(
          file.path(
            folder,
            "homerMotifs.motifs8"
          )
        ),
      
      homerMotifs_10 =
        file.exists(
          file.path(
            folder,
            "homerMotifs.motifs10"
          )
        ),
      
      homerMotifs_12 =
        file.exists(
          file.path(
            folder,
            "homerMotifs.motifs12"
          )
        ),
      
      motifFindingParameters =
        file.exists(
          file.path(
            folder,
            "motifFindingParameters.txt"
          )
        )
    ) %>%
      
      dplyr::bind_cols(
        meta
      )
  }
)


run_metadata <-
  drop_all_na_columns(
    run_metadata
  )


write.csv(
  run_metadata,
  file.path(
    analysis_root,
    "03_run_metadata.csv"
  ),
  row.names = FALSE
)


saveRDS(
  run_metadata,
  file.path(
    analysis_root,
    "03_run_metadata.rds"
  )
)


# ============================================================
# 6. READ KNOWN MOTIF RESULTS
# ============================================================

read_known_results <- function(folder) {
  
  f <- file.path(
    folder,
    "knownResults.txt"
  )
  
  if (!file.exists(f)) {
    return(tibble::tibble())
  }
  
  
  dat <- tryCatch(
    
    read.delim(
      f,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      fill = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    
    error = function(e) {
      
      warning(
        "Could not read: ",
        f,
        "\n",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  if (
    is.null(dat) ||
    nrow(dat) == 0
  ) {
    
    return(tibble::tibble())
    
  }
  
  
  # Preserve EXACT HOMER header
  original_headers <- names(dat)
  
  
  # Create safe names for programming
  names(dat) <-
    clean_colnames_unique(
      names(dat)
    )
  
  
  meta <- parse_homer_metadata(
    basename(folder)
  )
  
  
  dat <-
    dat %>%
    
    dplyr::mutate(
      
      result_type =
        "known",
      
      source_file =
        f,
      
      source_folder =
        folder,
      
      HOMER_original_headers =
        paste(
          original_headers,
          collapse = " | "
        )
      
    ) %>%
    
    dplyr::bind_cols(
      meta[
        rep(
          1,
          nrow(dat)
        ),
        ,
        drop = FALSE
      ]
    )
  
  
  dat
}


known_results <- purrr::map_dfr(
  homer_dirs,
  read_known_results
)


cat(
  "Known motif rows: ",
  nrow(known_results),
  "\n",
  sep = ""
)


write.csv(
  known_results,
  file.path(
    analysis_root,
    "01_known_motifs_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  known_results,
  file.path(
    analysis_root,
    "01_known_motifs_long.rds"
  )
)


# ============================================================
# 7. PRINT ACTUAL HOMER COLUMNS
# ============================================================

cat(
  "\n============================================\n"
)

cat(
  "HOMER columns detected\n"
)

cat(
  "============================================\n"
)

print(
  names(
    known_results
  )
)

cat(
  "============================================\n\n"
)


# ============================================================
# 8. COLUMN DETECTOR
# ============================================================

find_column <- function(
    data,
    patterns
) {
  
  nms <- names(data)
  
  hits <- which(
    
    vapply(
      nms,
      function(nm) {
        
        nm_lower <-
          tolower(nm)
        
        any(
          vapply(
            patterns,
            function(pattern) {
              
              grepl(
                pattern,
                nm_lower,
                perl = TRUE
              )
              
            },
            logical(1)
          )
        )
        
      },
      logical(1)
    )
  )
  
  
  if (
    length(hits) == 0
  ) {
    
    return(
      NA_character_
    )
    
  }
  
  
  nms[
    hits[1]
  ]
}


# ============================================================
# 9. DETECT HOMER COLUMNS
# ============================================================

p_col <- find_column(
  known_results,
  c(
    "^p_value$",
    "^pvalue$",
    "^p$",
    "p_value"
  )
)


q_col <- find_column(
  known_results,
  c(
    "q_value",
    "benjamini",
    "fdr"
  )
)


target_percent_col <- find_column(
  known_results,
  c(
    "target.*percent",
    "percent.*target",
    "target.*sequence.*percent"
  )
)


background_percent_col <- find_column(
  known_results,
  c(
    "background.*percent",
    "percent.*background",
    "background.*sequence.*percent"
  )
)


target_count_col <- find_column(
  known_results,
  c(
    "target.*sequence.*motif",
    "target.*motif",
    "target.*count"
  )
)


background_count_col <- find_column(
  known_results,
  c(
    "background.*sequence.*motif",
    "background.*motif",
    "background.*count"
  )
)


motif_col <- find_column(
  known_results,
  c(
    "^motif$",
    "motif_name",
    "motif"
  )
)


rank_col <- find_column(
  known_results,
  c(
    "^rank$"
  )
)


cat(
  "\n============================================\n"
)

cat(
  "Detected HOMER statistic columns\n"
)

cat(
  "============================================\n"
)

cat(
  "P-value:          ",
  p_col,
  "\n",
  sep = ""
)

cat(
  "Q-value:          ",
  q_col,
  "\n",
  sep = ""
)

cat(
  "Target %:         ",
  target_percent_col,
  "\n",
  sep = ""
)

cat(
  "Background %:     ",
  background_percent_col,
  "\n",
  sep = ""
)

cat(
  "Target count:     ",
  target_count_col,
  "\n",
  sep = ""
)

cat(
  "Background count: ",
  background_count_col,
  "\n",
  sep = ""
)

cat(
  "Motif:            ",
  motif_col,
  "\n",
  sep = ""
)

cat(
  "Rank:             ",
  rank_col,
  "\n",
  sep = ""
)

cat(
  "============================================\n\n"
)


# ============================================================
# 10. STANDARDIZE HOMER STATISTICS
# ============================================================

known_results_std <-
  known_results %>%
  
  dplyr::mutate(
    
    homer_pvalue =
      if (!is.na(p_col)) {
        
        safe_numeric(
          .data[[p_col]]
        )
        
      } else {
        
        NA_real_
        
      },
    
    
    homer_qvalue =
      if (!is.na(q_col)) {
        
        safe_numeric(
          .data[[q_col]]
        )
        
      } else {
        
        NA_real_
        
      },
    
    
    homer_target_percent =
      if (!is.na(target_percent_col)) {
        
        safe_numeric(
          .data[[target_percent_col]]
        )
        
      } else {
        
        NA_real_
        
      },
    
    
    homer_background_percent =
      if (!is.na(background_percent_col)) {
        
        safe_numeric(
          .data[[background_percent_col]]
        )
        
      } else {
        
        NA_real_
        
      },
    
    
    homer_target_count =
      if (!is.na(target_count_col)) {
        
        safe_numeric(
          .data[[target_count_col]]
        )
        
      } else {
        
        NA_real_
        
      },
    
    
    homer_background_count =
      if (!is.na(background_count_col)) {
        
        safe_numeric(
          .data[[background_count_col]]
        )
        
      } else {
        
        NA_real_
        
      }
    
  )


# ============================================================
# 11. DERIVED TARGET/BACKGROUND RATIO
# ============================================================
#
# This is NOT claimed to be a native HOMER statistic.
#
# Preferred:
#   target percentage / background percentage
#
# Fallback:
#   target count / background count
#
# ============================================================

known_results_std <-
  known_results_std %>%
  
  dplyr::mutate(
    
    target_background_ratio =
      
      dplyr::case_when(
        
        is.finite(
          homer_target_percent
        ) &
          
          is.finite(
            homer_background_percent
          ) &
          
          homer_background_percent > 0
        
        ~
          
          homer_target_percent /
          homer_background_percent,
        
        
        is.finite(
          homer_target_count
        ) &
          
          is.finite(
            homer_background_count
          ) &
          
          homer_background_count > 0
        
        ~
          
          homer_target_count /
          homer_background_count,
        
        
        TRUE ~
          
          NA_real_
        
      ),
    
    
    log2_target_background =
      
      dplyr::case_when(
        
        is.finite(
          target_background_ratio
        ) &
          
          target_background_ratio > 0
        
        ~
          
          log2(
            target_background_ratio
          ),
        
        
        TRUE ~
          
          NA_real_
        
      ),
    
    
    minus_log10_p =
      
      dplyr::case_when(
        
        is.finite(
          homer_pvalue
        ) &
          
          homer_pvalue > 0
        
        ~
          
          -log10(
            homer_pvalue
          ),
        
        
        homer_pvalue == 0
        
        ~
          
          300,
        
        
        TRUE ~
          
          NA_real_
        
      ),
    
    
    significant_FDR =
      
      !is.na(
        homer_qvalue
      ) &
      
      is.finite(
        homer_qvalue
      ) &
      
      homer_qvalue <
      FDR_THRESHOLD,
    
    
    strong_enrichment =
      
      !is.na(
        homer_qvalue
      ) &
      
      is.finite(
        homer_qvalue
      ) &
      
      homer_qvalue <
      FDR_THRESHOLD &
      
      is.finite(
        target_background_ratio
      ) &
      
      target_background_ratio >=
      MIN_TARGET_BACKGROUND_RATIO
    
  )


# Remove 100% NA columns
known_results_std <-
  drop_all_na_columns(
    known_results_std
  )


write.csv(
  known_results_std,
  file.path(
    analysis_root,
    "01_known_motifs_standardized_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  known_results_std,
  file.path(
    analysis_root,
    "01_known_motifs_standardized_long.rds"
  )
)


# ============================================================
# 12. READ DE NOVO HOMER RESULTS
# ============================================================

read_denovo_results <- function(folder) {
  
  f <- file.path(
    folder,
    "homerResults.html"
  )
  
  
  if (!file.exists(f)) {
    return(
      tibble::tibble()
    )
  }
  
  
  page <- tryCatch(
    
    rvest::read_html(
      f
    ),
    
    error = function(e) {
      
      warning(
        "Could not read ",
        f,
        "\n",
        conditionMessage(e)
      )
      
      NULL
      
    }
  )
  
  
  if (is.null(page)) {
    return(
      tibble::tibble()
    )
  }
  
  
  tables <- tryCatch(
    
    rvest::html_table(
      page,
      fill = TRUE
    ),
    
    error = function(e) {
      
      list()
      
    }
  )
  
  
  if (
    length(tables) == 0
  ) {
    
    warning(
      "No HTML tables found in ",
      f
    )
    
    return(
      tibble::tibble()
    )
  }
  
  
  # Score tables to identify likely de novo table
  scores <- purrr::map_int(
    
    tables,
    
    function(tb) {
      
      nm <- tolower(
        names(tb)
      )
      
      score <- 0
      
      
      if (
        any(
          grepl(
            "rank",
            nm
          )
        )
      ) {
        
        score <- score + 1
        
      }
      
      
      if (
        any(
          grepl(
            "p.?value|pvalue",
            nm,
            perl = TRUE
          )
        )
      ) {
        
        score <- score + 1
        
      }
      
      
      if (
        any(
          grepl(
            "target",
            nm
          )
        )
      ) {
        
        score <- score + 1
        
      }
      
      
      if (
        any(
          grepl(
            "background",
            nm
          )
        )
      ) {
        
        score <- score + 1
        
      }
      
      
      if (
        any(
          grepl(
            "known|match|similar|factor",
            nm
          )
        )
      ) {
        
        score <- score + 1
        
      }
      
      
      score
      
    }
  )
  
  
  best <- which.max(
    scores
  )
  
  
  if (
    length(best) == 0 ||
    scores[best] < 2
  ) {
    
    warning(
      "Could not identify a likely de novo motif table in ",
      f
    )
    
    return(
      tibble::tibble()
    )
  }
  
  
  dat <- tables[[best]]
  
  
  names(dat) <-
    clean_colnames_unique(
      names(dat)
    )
  
  
  meta <-
    parse_homer_metadata(
      basename(folder)
    )
  
  
  dat <-
    dat %>%
    
    dplyr::mutate(
      
      result_type =
        "de_novo",
      
      source_file =
        f,
      
      source_folder =
        folder
      
    ) %>%
    
    dplyr::bind_cols(
      meta[
        rep(
          1,
          nrow(dat)
        ),
        ,
        drop = FALSE
      ]
    )
  
  
  dat
}


denovo_results <-
  purrr::map_dfr(
    homer_dirs,
    read_denovo_results
  )


denovo_results <-
  drop_all_na_columns(
    denovo_results
  )


cat(
  "De novo motif rows: ",
  nrow(denovo_results),
  "\n",
  sep = ""
)


write.csv(
  denovo_results,
  file.path(
    analysis_root,
    "02_denovo_motifs_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  denovo_results,
  file.path(
    analysis_root,
    "02_denovo_motifs_long.rds"
  )
)


# ============================================================
# 13. READ HOMER PARAMETERS
# ============================================================

read_parameters <- function(folder) {
  
  f <- file.path(
    folder,
    "motifFindingParameters.txt"
  )
  
  
  if (!file.exists(f)) {
    
    return(
      tibble::tibble()
    )
    
  }
  
  
  txt <- readLines(
    f,
    warn = FALSE
  )
  
  
  if (
    length(txt) == 0
  ) {
    
    return(
      tibble::tibble()
    )
  }
  
  
  meta <-
    parse_homer_metadata(
      basename(folder)
    )
  
  
  tibble::tibble(
    
    parameter_line =
      seq_along(txt),
    
    parameter_text =
      txt,
    
    source_folder =
      folder
    
  ) %>%
    
    dplyr::bind_cols(
      meta[
        rep(
          1,
          length(txt)
        ),
        ,
        drop = FALSE
      ]
    )
}


parameters_long <-
  purrr::map_dfr(
    homer_dirs,
    read_parameters
  )


parameters_long <-
  drop_all_na_columns(
    parameters_long
  )


write.csv(
  parameters_long,
  file.path(
    analysis_root,
    "04_parameters_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  parameters_long,
  file.path(
    analysis_root,
    "04_parameters_long.rds"
  )
)


# ============================================================
# 14. FINAL MOTIF KEY
# ============================================================

if (!is.na(motif_col)) {
  
  known_results_std <-
    known_results_std %>%
    
    dplyr::mutate(
      
      motif_key =
        as.character(
          .data[[motif_col]]
        )
      
    )
  
} else {
  
  known_results_std$motif_key <-
    paste0(
      "motif_",
      seq_len(
        nrow(
          known_results_std
        )
      )
    )
}


# Save again because motif_key was added
known_results_std <-
  drop_all_na_columns(
    known_results_std
  )


write.csv(
  known_results_std,
  file.path(
    analysis_root,
    "01_known_motifs_standardized_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  known_results_std,
  file.path(
    analysis_root,
    "01_known_motifs_standardized_long.rds"
  )
)


# ============================================================
# 15. SUMMARY BY HOMER RUN
# ============================================================

run_summary <-
  known_results_std %>%
  
  dplyr::group_by(
    dataset,
    analysis,
    cell_type,
    stage,
    run_name
  ) %>%
  
  dplyr::summarise(
    
    n_known_motifs =
      dplyr::n(),
    
    n_FDR05 =
      sum(
        significant_FDR,
        na.rm = TRUE
      ),
    
    n_strong_enrichment =
      sum(
        strong_enrichment,
        na.rm = TRUE
      ),
    
    best_qvalue =
      safe_min(
        homer_qvalue
      ),
    
    strongest_target_background_ratio =
      safe_max(
        target_background_ratio
      ),
    
    .groups = "drop"
  )


run_summary <-
  drop_all_na_columns(
    run_summary
  )


write.csv(
  run_summary,
  file.path(
    analysis_root,
    "10_summary_by_run.csv"
  ),
  row.names = FALSE
)


saveRDS(
  run_summary,
  file.path(
    analysis_root,
    "10_summary_by_run.rds"
  )
)


# ============================================================
# 16. MOTIF RECURRENCE ACROSS ALL RUNS
# ============================================================

motif_recurrence <-
  known_results_std %>%
  
  dplyr::group_by(
    motif_key
  ) %>%
  
  dplyr::summarise(
    
    motif_name =
      dplyr::first(
        motif_key
      ),
    
    n_runs =
      dplyr::n_distinct(
        run_name
      ),
    
    n_significant_runs =
      dplyr::n_distinct(
        run_name[
          significant_FDR
        ]
      ),
    
    n_strong_runs =
      dplyr::n_distinct(
        run_name[
          strong_enrichment
        ]
      ),
    
    n_datasets =
      dplyr::n_distinct(
        dataset[
          significant_FDR
        ]
      ),
    
    significant_datasets =
      paste(
        sort(
          unique(
            dataset[
              significant_FDR &
                !is.na(dataset)
            ]
          )
        ),
        collapse = ";"
      ),
    
    significant_cell_types =
      paste(
        sort(
          unique(
            cell_type[
              significant_FDR &
                !is.na(cell_type)
            ]
          )
        ),
        collapse = ";"
      ),
    
    significant_stages =
      paste(
        sort(
          unique(
            stage[
              significant_FDR &
                !is.na(stage)
            ]
          )
        ),
        collapse = ";"
      ),
    
    median_target_background_ratio =
      safe_median(
        target_background_ratio[
          significant_FDR
        ]
      ),
    
    max_target_background_ratio =
      safe_max(
        target_background_ratio[
          significant_FDR
        ]
      ),
    
    best_qvalue =
      safe_min(
        homer_qvalue
      ),
    
    .groups = "drop"
  ) %>%
  
  dplyr::arrange(
    dplyr::desc(
      n_strong_runs
    ),
    dplyr::desc(
      n_significant_runs
    ),
    best_qvalue
  )


motif_recurrence <-
  drop_all_na_columns(
    motif_recurrence
  )


write.csv(
  motif_recurrence,
  file.path(
    analysis_root,
    "05_motif_recurrence.csv"
  ),
  row.names = FALSE
)


saveRDS(
  motif_recurrence,
  file.path(
    analysis_root,
    "05_motif_recurrence.rds"
  )
)


# ============================================================
# 17. PRENATAL VS POSTNATAL COMPARISON
# ============================================================
#
# This is DESCRIPTIVE.
#
# It does NOT claim that significance differs statistically
# between prenatal and postnatal.
# ============================================================

stage_known <-
  known_results_std %>%
  
  dplyr::filter(
    stage %in% c(
      "prenatal",
      "postnatal"
    )
  ) %>%
  
  dplyr::group_by(
    dataset,
    analysis,
    cell_type,
    motif_key
  ) %>%
  
  dplyr::summarise(
    
    prenatal_ratio = {
      
      x <-
        target_background_ratio[
          stage == "prenatal"
        ]
      
      safe_max(x)
      
    },
    
    postnatal_ratio = {
      
      x <-
        target_background_ratio[
          stage == "postnatal"
        ]
      
      safe_max(x)
      
    },
    
    prenatal_q = {
      
      x <-
        homer_qvalue[
          stage == "prenatal"
        ]
      
      safe_min(x)
      
    },
    
    postnatal_q = {
      
      x <-
        homer_qvalue[
          stage == "postnatal"
        ]
      
      safe_min(x)
      
    },
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    log2_prenatal_postnatal =
      dplyr::case_when(
        
        is.finite(
          prenatal_ratio
        ) &
          
          is.finite(
            postnatal_ratio
          ) &
          
          prenatal_ratio > 0 &
          
          postnatal_ratio > 0
        
        ~
          
          log2(
            prenatal_ratio /
              postnatal_ratio
          ),
        
        TRUE ~
          
          NA_real_
        
      ),
    
    
    stage_pattern =
      dplyr::case_when(
        
        prenatal_q < FDR_THRESHOLD &
          
          postnatal_q < FDR_THRESHOLD
        
        ~
          
          "significant_both",
        
        
        prenatal_q < FDR_THRESHOLD &
          
          (
            is.na(postnatal_q) |
              postnatal_q >= FDR_THRESHOLD
          )
        
        ~
          
          "prenatal_only_significant",
        
        
        postnatal_q < FDR_THRESHOLD &
          
          (
            is.na(prenatal_q) |
              prenatal_q >= FDR_THRESHOLD
          )
        
        ~
          
          "postnatal_only_significant",
        
        
        is.finite(
          log2_prenatal_postnatal
        ) &
          
          log2_prenatal_postnatal > 0
        
        ~
          
          "prenatal_higher_ratio",
        
        
        is.finite(
          log2_prenatal_postnatal
        ) &
          
          log2_prenatal_postnatal < 0
        
        ~
          
          "postnatal_higher_ratio",
        
        
        TRUE ~
          
          "insufficient_data"
        
      )
  )


stage_known <-
  drop_all_na_columns(
    stage_known
  )


write.csv(
  stage_known,
  file.path(
    analysis_root,
    "06_motif_stage_shift.csv"
  ),
  row.names = FALSE
)


saveRDS(
  stage_known,
  file.path(
    analysis_root,
    "06_motif_stage_shift.rds"
  )
)


# ============================================================
# 18. CELL-TYPE SPECIFICITY
# ============================================================

celltype_specificity <-
  known_results_std %>%
  
  dplyr::filter(
    significant_FDR
  ) %>%
  
  dplyr::group_by(
    dataset,
    analysis,
    motif_key
  ) %>%
  
  dplyr::summarise(
    
    n_cell_types =
      dplyr::n_distinct(
        cell_type[
          !is.na(cell_type)
        ]
      ),
    
    
    cell_types =
      paste(
        sort(
          unique(
            cell_type[
              !is.na(cell_type)
            ]
          )
        ),
        collapse = ";"
      ),
    
    
    best_cell_type = {
      
      valid_idx <-
        which(
          !is.na(
            target_background_ratio
          ) &
            
            is.finite(
              target_background_ratio
            ) &
            
            !is.na(
              cell_type
            )
        )
      
      
      if (
        length(valid_idx) == 0
      ) {
        
        NA_character_
        
      } else {
        
        cell_type[
          valid_idx[
            which.max(
              target_background_ratio[
                valid_idx
              ]
            )
          ]
        ]
        
      }
      
    },
    
    
    best_target_background_ratio = {
      
      x <-
        target_background_ratio[
          is.finite(
            target_background_ratio
          )
        ]
      
      safe_max(x)
      
    },
    
    
    median_target_background_ratio = {
      
      x <-
        target_background_ratio[
          is.finite(
            target_background_ratio
          )
        ]
      
      safe_median(x)
      
    },
    
    
    .groups = "drop"
  ) %>%
  
  dplyr::mutate(
    
    specificity_ratio =
      
      dplyr::case_when(
        
        is.finite(
          best_target_background_ratio
        ) &
          
          is.finite(
            median_target_background_ratio
          ) &
          
          median_target_background_ratio > 0
        
        ~
          
          best_target_background_ratio /
          median_target_background_ratio,
        
        TRUE ~
          
          NA_real_
        
      )
    
  ) %>%
  
  dplyr::arrange(
    dataset,
    dplyr::desc(
      specificity_ratio
    )
  )


celltype_specificity <-
  drop_all_na_columns(
    celltype_specificity
  )


write.csv(
  celltype_specificity,
  file.path(
    analysis_root,
    "07_celltype_specificity.csv"
  ),
  row.names = FALSE
)


saveRDS(
  celltype_specificity,
  file.path(
    analysis_root,
    "07_celltype_specificity.rds"
  )
)


# ============================================================
# 19. DE NOVO SUMMARY
# ============================================================

if (
  nrow(denovo_results) > 0
) {
  
  denovo_summary <-
    denovo_results %>%
    
    dplyr::group_by(
      dataset,
      analysis,
      cell_type,
      stage,
      run_name
    ) %>%
    
    dplyr::summarise(
      
      n_denovo_motifs =
        dplyr::n(),
      
      .groups = "drop"
      
    )
  
  
  denovo_summary <-
    drop_all_na_columns(
      denovo_summary
    )
  
  
  write.csv(
    denovo_summary,
    file.path(
      analysis_root,
      "08_denovo_summary.csv"
    ),
    row.names = FALSE
  )
  
  
  saveRDS(
    denovo_summary,
    file.path(
      analysis_root,
      "08_denovo_summary.rds"
    )
  )
  
} else {
  
  denovo_summary <-
    tibble::tibble(
      message =
        "No de novo result tables were successfully extracted."
    )
}


# ============================================================
# 20. DE NOVO -> KNOWN CONVERGENCE
# ============================================================

if (
  nrow(denovo_results) > 0
) {
  
  dn_names <-
    names(
      denovo_results
    )
  
  
  match_candidates <-
    dn_names[
      vapply(
        dn_names,
        function(x) {
          
          grepl(
            "known|match|similar|factor",
            tolower(x)
          )
          
        },
        logical(1)
      )
    ]
  
  
  if (
    length(match_candidates) > 0
  ) {
    
    best_match_col <-
      match_candidates[1]
    
    
    denovo_known_convergence <-
      denovo_results %>%
      
      dplyr::mutate(
        
        best_known_match =
          as.character(
            .data[[best_match_col]]
          )
        
      ) %>%
      
      dplyr::filter(
        
        !is.na(
          best_known_match
        ),
        
        nzchar(
          best_known_match
        )
        
      )
    
    
  } else {
    
    denovo_known_convergence <-
      tibble::tibble()
    
  }
  
  
} else {
  
  denovo_known_convergence <-
    tibble::tibble()
  
}


denovo_known_convergence <-
  drop_all_na_columns(
    denovo_known_convergence
  )


write.csv(
  denovo_known_convergence,
  file.path(
    analysis_root,
    "09_denovo_known_convergence.csv"
  ),
  row.names = FALSE
)


saveRDS(
  denovo_known_convergence,
  file.path(
    analysis_root,
    "09_denovo_known_convergence.rds"
  )
)


# ============================================================
# 21. MOTIF HEATMAP
# ============================================================
#
# Uses ONLY finite derived log2 ratios.
#
# No NA / Inf / -Inf values are sent into pheatmap.
# ============================================================

heatmap_long <-
  known_results_std %>%
  
  dplyr::filter(
    
    is.finite(
      log2_target_background
    )
    
  ) %>%
  
  dplyr::group_by(
    run_name,
    motif_key
  ) %>%
  
  dplyr::summarise(
    
    value =
      median(
        log2_target_background,
        na.rm = TRUE
      ),
    
    .groups = "drop"
    
  )


if (
  nrow(heatmap_long) > 0
) {
  
  heatmap_wide <-
    heatmap_long %>%
    
    tidyr::pivot_wider(
      
      names_from =
        motif_key,
      
      values_from =
        value
      
    )
  
  
  heatmap_wide <-
    drop_all_na_columns(
      heatmap_wide
    )
  
  
  write.csv(
    heatmap_wide,
    file.path(
      analysis_root,
      "11_motif_heatmap_matrix.csv"
    ),
    row.names = FALSE
  )
  
  
  if (
    nrow(heatmap_wide) >= 2 &&
    ncol(heatmap_wide) >= 3
  ) {
    
    heat_mat <-
      as.matrix(
        heatmap_wide[
          ,
          -1,
          drop = FALSE
        ]
      )
    
    
    rownames(heat_mat) <-
      heatmap_wide$run_name
    
    
    heat_mat[
      !is.finite(
        heat_mat
      )
    ] <- NA_real_
    
    
    keep_rows <-
      apply(
        heat_mat,
        1,
        function(x) {
          
          any(
            is.finite(x)
          )
          
        }
      )
    
    
    keep_cols <-
      apply(
        heat_mat,
        2,
        function(x) {
          
          any(
            is.finite(x)
          )
          
        }
      )
    
    
    heat_mat <-
      heat_mat[
        keep_rows,
        keep_cols,
        drop = FALSE
      ]
    
    
    if (
      nrow(heat_mat) >= 2 &&
      ncol(heat_mat) >= 2
    ) {
      
      pdf(
        file.path(
          analysis_root,
          "11_motif_enrichment_heatmap.pdf"
        ),
        width = 16,
        height = 12
      )
      
      
      pheatmap::pheatmap(
        
        heat_mat,
        
        cluster_rows =
          TRUE,
        
        cluster_cols =
          TRUE,
        
        scale =
          "column",
        
        na_col =
          "grey90",
        
        fontsize =
          7,
        
        main =
          "HOMER motif target/background ratio"
        
      )
      
      
      dev.off()
      
      
    } else {
      
      message(
        "Heatmap skipped: insufficient usable rows/columns."
      )
      
    }
    
    
  } else {
    
    message(
      "Heatmap skipped: matrix too small."
    )
    
  }
  
  
} else {
  
  heatmap_wide <-
    tibble::tibble(
      
      message =
        "No finite derived enrichment values were available."
      
    )
  
  
  write.csv(
    heatmap_wide,
    file.path(
      analysis_root,
      "11_motif_heatmap_matrix.csv"
    ),
    row.names = FALSE
  )
  
  
  message(
    "Heatmap skipped: no finite derived enrichment values."
  )
  
}


# ============================================================
# 22. MOTIF RANK STABILITY
# ============================================================

if (
  !is.na(q_col)
) {
  
  rank_stability <-
    known_results_std %>%
    
    dplyr::group_by(
      dataset,
      analysis,
      cell_type,
      stage
    ) %>%
    
    dplyr::arrange(
      homer_qvalue,
      .by_group = TRUE
    ) %>%
    
    dplyr::mutate(
      
      within_run_rank =
        dplyr::row_number()
      
    ) %>%
    
    dplyr::ungroup()
  
  
} else {
  
  rank_stability <-
    tibble::tibble()
  
}


rank_stability <-
  drop_all_na_columns(
    rank_stability
  )


write.csv(
  rank_stability,
  file.path(
    analysis_root,
    "12_motif_rank_stability_long.csv"
  ),
  row.names = FALSE
)


saveRDS(
  rank_stability,
  file.path(
    analysis_root,
    "12_motif_rank_stability_long.rds"
  )
)


# ============================================================
# 23. OPTIONAL MOTIF INSTANCE ANALYSIS
# ============================================================
#
# This is OFF by default.
#
# Once the main results are validated, set:
#
# RUN_MOTIF_INSTANCES <- TRUE
#
# ============================================================

run_one_instance_analysis <- function(folder) {
  
  base_name <-
    stringr::str_remove(
      basename(folder),
      "_HOMER$"
    )
  
  
  peak_files <-
    list.files(
      homer_root,
      pattern =
        paste0(
          "^",
          stringr::str_escape(
            base_name
          ),
          "\\.homerpeaks\\.txt$"
        ),
      full.names = TRUE
    )
  
  
  if (
    length(peak_files) == 0
  ) {
    
    warning(
      "No .homerpeaks.txt file found for ",
      base_name
    )
    
    return(NULL)
  }
  
  
  peak_file <-
    peak_files[1]
  
  
  motif_dir <-
    file.path(
      folder,
      "homerResults"
    )
  
  
  motif_files <-
    list.files(
      motif_dir,
      pattern =
        "\\.motif$",
      full.names = TRUE
    )
  
  
  if (
    length(motif_files) == 0
  ) {
    
    warning(
      "No de novo .motif files found in ",
      motif_dir
    )
    
    return(NULL)
  }
  
  
  # First 20 de novo motifs
  motif_files <-
    head(
      motif_files,
      20
    )
  
  
  output_file <-
    file.path(
      instance_dir,
      paste0(
        base_name,
        "_motif_instances.txt"
      )
    )
  
  
  args <- c(
    
    peak_file,
    
    HOMER_GENOME,
    
    "-m",
    
    motif_files
    
  )
  
  
  ok <- tryCatch(
    
    {
      
      result <-
        system2(
          HOMER_ANNOTATE,
          args = args,
          stdout = output_file,
          stderr = TRUE
        )
      
      isTRUE(
        result == 0
      )
      
    },
    
    error = function(e) {
      
      warning(
        "annotatePeaks.pl failed for ",
        base_name,
        ": ",
        conditionMessage(e)
      )
      
      FALSE
    }
  )
  
  
  if (ok) {
    
    return(
      output_file
    )
    
  }
  
  
  NULL
}


if (
  RUN_MOTIF_INSTANCES
) {
  
  message(
    "\nRunning motif-instance analysis...\n"
  )
  
  
  instance_files <-
    purrr::map(
      homer_dirs,
      run_one_instance_analysis
    )
  
  
  instance_files <-
    unlist(
      instance_files
    )
  
  
  writeLines(
    instance_files,
    file.path(
      analysis_root,
      "motif_instance_files.txt"
    )
  )
  
}


# ============================================================
# 24. SAVE EXCEL WORKBOOK
# ============================================================

wb <-
  openxlsx::createWorkbook()


add_data_sheet <- function(
    wb,
    sheet_name,
    data
) {
  
  openxlsx::addWorksheet(
    wb,
    sheet_name
  )
  
  
  if (
    is.null(data) ||
    nrow(data) == 0
  ) {
    
    openxlsx::writeData(
      wb,
      sheet_name,
      data.frame(
        message =
          "No data available"
      )
    )
    
  } else {
    
    openxlsx::writeData(
      wb,
      sheet_name,
      data,
      withFilter = TRUE
    )
  }
  
}


# Main tables

add_data_sheet(
  wb,
  "KNOWN_LONG",
  known_results_std
)


add_data_sheet(
  wb,
  "DENOVO_LONG",
  denovo_results
)


add_data_sheet(
  wb,
  "RUN_SUMMARY",
  run_summary
)


add_data_sheet(
  wb,
  "MOTIF_RECURRENCE",
  motif_recurrence
)


add_data_sheet(
  wb,
  "STAGE_SHIFT",
  stage_known
)


add_data_sheet(
  wb,
  "CELLTYPE_SPECIFICITY",
  celltype_specificity
)


add_data_sheet(
  wb,
  "DENOVO_SUMMARY",
  denovo_summary
)


add_data_sheet(
  wb,
  "DENOVO_CONVERGENCE",
  denovo_known_convergence
)


add_data_sheet(
  wb,
  "PARAMETERS",
  parameters_long
)


add_data_sheet(
  wb,
  "RUN_METADATA",
  run_metadata
)


add_data_sheet(
  wb,
  "HEATMAP_MATRIX",
  heatmap_wide
)


add_data_sheet(
  wb,
  "RANK_STABILITY",
  rank_stability
)


# ============================================================
# 25. EXCEL FORMATTING
# ============================================================

header_style <-
  openxlsx::createStyle(
    
    textDecoration =
      "bold",
    
    fgFill =
      "#1F4E78",
    
    fontColour =
      "#FFFFFF",
    
    halign =
      "center",
    
    valign =
      "center"
    
  )


for (
  sheet in names(wb)
) {
  
  dat <-
    openxlsx::readWorkbook(
      wb,
      sheet
    )
  
  
  if (
    ncol(dat) > 0
  ) {
    
    openxlsx::addStyle(
      
      wb,
      
      sheet,
      
      header_style,
      
      rows = 1,
      
      cols =
        seq_len(
          ncol(dat)
        ),
      
      gridExpand =
        TRUE
      
    )
    
    
    openxlsx::freezePane(
      
      wb,
      
      sheet,
      
      firstRow =
        TRUE
      
    )
    
    
    openxlsx::setColWidths(
      
      wb,
      
      sheet,
      
      cols =
        seq_len(
          ncol(dat)
        ),
      
      widths =
        "auto"
      
    )
  }
}


# Scientific notation
# for HOMER P/Q values

for (
  sheet in c(
    "KNOWN_LONG",
    "RUN_SUMMARY",
    "MOTIF_RECURRENCE"
  )
) {
  
  if (
    sheet %in%
    names(wb)
  ) {
    
    dat <-
      openxlsx::readWorkbook(
        wb,
        sheet
      )
    
    
    for (
      col_name in c(
        "homer_pvalue",
        "homer_qvalue",
        "best_qvalue"
      )
    ) {
      
      idx <-
        which(
          names(dat) ==
            col_name
        )
      
      
      if (
        length(idx) == 1 &&
        nrow(dat) > 0
      ) {
        
        openxlsx::addStyle(
          
          wb,
          
          sheet,
          
          openxlsx::createStyle(
            numFmt =
              "0.00E+00"
          ),
          
          rows =
            2:(nrow(dat) + 1),
          
          cols =
            idx,
          
          gridExpand =
            TRUE
          
        )
        
      }
      
    }
  }
}


# ============================================================
# 26. SAVE EXCEL
# ============================================================

xlsx_out <-
  file.path(
    analysis_root,
    "ALL_HOMER_RESULTS.xlsx"
  )


openxlsx::saveWorkbook(
  
  wb,
  
  xlsx_out,
  
  overwrite =
    TRUE
  
)


# ============================================================
# 27. FINAL VALIDATION
# ============================================================

cat("\n")
cat("====================================================\n")
cat("HOMER INTEGRATED ANALYSIS COMPLETE\n")
cat("====================================================\n")

cat(
  "HOMER result folders: ",
  length(homer_dirs),
  "\n",
  sep = ""
)

cat(
  "Known motif rows: ",
  nrow(known_results_std),
  "\n",
  sep = ""
)

cat(
  "De novo rows: ",
  nrow(denovo_results),
  "\n",
  sep = ""
)

cat(
  "Significant known-motif rows: ",
  sum(
    known_results_std$significant_FDR,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "\nOutput folder:\n",
  analysis_root,
  "\n",
  sep = ""
)

cat(
  "\nExcel workbook:\n",
  xlsx_out,
  "\n",
  sep = ""
)

cat(
  "\n====================================================\n"
)

# ============================================================
# END
# ============================================================