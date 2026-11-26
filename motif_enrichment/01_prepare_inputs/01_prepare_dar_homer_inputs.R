# ============================================================
# Generate top-N peaks per analysis + HOMER-format input files
# Fully self-contained: reads the Excel sheets, builds the top-N peak
# sets, and writes HOMER peak files + a run script. No external state
# assumed.
# ============================================================

library(GenomicRanges)
library(IRanges)
library(dplyr)
library(readxl)
library(conflicted)

conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::slice)
conflicts_prefer(dplyr::desc)
conflicts_prefer(GenomicRanges::reduce)

# ---- CONFIG ----
top_n_peaks <- 2000

# Source Excel files live in Dropbox -- fine, we only READ from here
input_dir <- "/Users/nrohani/Library/CloudStorage/Dropbox-Izzy.barraclough/Narjes Rohani/BrainVar_NRohani/Supplementary_Data"

# Everything WRITTEN (peak files, run script, and later HOMER's own output
# folders + temp files) goes to a local, non-synced directory. Dropbox's
# sync daemon races with HOMER's fast-churning temp files and causes
# "file not found" errors mid-run if outputs live inside a synced folder.
out_dir <- path.expand("~/homer_runs/motif_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sheet_rules <- list(
  list(workbook = "Data S1.xlsx", sheet = "DAR_broad_cell_type_prenatal",  group_cols = c("cluster")),
  list(workbook = "Data S1.xlsx", sheet = "DAR_broad_cell_type_postnatal",  group_cols = c("cluster")),
  list(workbook = "Data S1.xlsx", sheet = "DAR_neuron_cell_type_prenatal",  group_cols = c("cluster")),
  list(workbook = "Data S1.xlsx", sheet = "DAR_neuron_cell_type_postnatal", group_cols = c("cluster")),
  list(workbook = "Data S2.xlsx", sheet = "DAR_development",               group_cols = c("cell_type", "stage")),
  list(workbook = "Data S3.xlsx", sheet = "DAR_celltype&development",      group_cols = c("cell_type", "stage"))
)

# ============================================================
# ---- helpers: parsing / cleaning ----
# ============================================================

clean_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

normalize_chr <- function(chr) {
  chr <- as.character(chr)
  chr <- ifelse(is.na(chr), NA_character_, chr)
  chr <- ifelse(startsWith(tolower(chr), "chr"), chr, paste0("chr", chr))
  chr
}

parse_peak_vector <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (length(x) == 0) return(GRanges())
  
  toks <- unlist(strsplit(x, "[;|,\\n\\r]+", perl = TRUE))
  toks <- trimws(toks)
  toks <- toks[nzchar(toks)]
  if (length(toks) == 0) return(GRanges())
  
  gr_list <- vector("list", length(toks))
  keep <- logical(length(toks))
  
  for (i in seq_along(toks)) {
    tok <- toks[i]
    
    m1 <- regexec("^(.+?):([0-9]+)-([0-9]+)$", tok)
    r1 <- regmatches(tok, m1)[[1]]
    if (length(r1) == 4) {
      chr <- normalize_chr(r1[2])
      gr_list[[i]] <- GRanges(seqnames = chr, ranges = IRanges(start = as.integer(r1[3]), end = as.integer(r1[4])))
      keep[i] <- TRUE
      next
    }
    
    parts <- unlist(strsplit(tok, "[[:space:]]+", perl = TRUE))
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 3 && grepl("^[0-9]+$", parts[2]) && grepl("^[0-9]+$", parts[3])) {
      chr <- normalize_chr(parts[1])
      gr_list[[i]] <- GRanges(seqnames = chr, ranges = IRanges(start = as.integer(parts[2]), end = as.integer(parts[3])))
      keep[i] <- TRUE
      next
    }
  }
  
  if (!any(keep)) return(GRanges())
  do.call(c, gr_list[keep])
}

pick_rank_cols <- function(df) {
  nm_lower <- tolower(names(df))
  
  fc_idx <- which(
    grepl("log2.*fold|fold.*change|foldchange|\\bfc\\b", nm_lower) |
      nm_lower %in% c("log2foldchange", "foldchange", "fold_change", "log2fc", "fc")
  )
  
  padj_idx <- which(
    grepl("padj|adj.*p|qval|fdr", nm_lower) |
      nm_lower %in% c("padj", "adj.p.val", "adj_pval", "qval", "fdr", "adj_p_value")
  )
  
  fc_col   <- if (length(fc_idx) > 0) names(df)[fc_idx[1]] else NA_character_
  padj_col <- if (length(padj_idx) > 0) names(df)[padj_idx[1]] else NA_character_
  
  list(fc = fc_col, padj = padj_col)
}

extract_sets_from_sheet <- function(file, sheet, group_cols = NULL, top_n_peaks = 2000) {
  df <- read_xlsx(file, sheet = sheet, col_types = "text")
  if (nrow(df) == 0) return(list())
  
  nm_lower <- tolower(names(df))
  coord_idx <- match(c("seqnames", "start", "end"), nm_lower)
  
  if (all(!is.na(coord_idx)) && !is.null(group_cols)) {
    group_idx <- match(tolower(group_cols), nm_lower)
    if (all(!is.na(group_idx))) {
      coord_cols <- names(df)[coord_idx]
      group_cols_actual <- names(df)[group_idx]
      
      score_cols <- pick_rank_cols(df)
      
      keep_cols <- unique(c(coord_cols, group_cols_actual, score_cols$fc, score_cols$padj))
      keep_cols <- keep_cols[!is.na(keep_cols)]
      df2 <- df[, keep_cols, drop = FALSE]
      df2 <- df2[complete.cases(df2[, coord_cols]), , drop = FALSE]
      
      # rank and keep only top N peaks per group if score columns exist
      if (!is.na(score_cols$fc) && !is.na(score_cols$padj)) {
        message(
          "Ranking peaks in ", sheet,
          " using ", score_cols$fc, " (desc) and ", score_cols$padj, " (asc); keeping top ",
          top_n_peaks, " per group"
        )
        
        df2[[score_cols$fc]]   <- suppressWarnings(as.numeric(df2[[score_cols$fc]]))
        df2[[score_cols$padj]] <- suppressWarnings(as.numeric(df2[[score_cols$padj]]))
        
        df2 <- df2 %>%
          dplyr::group_by(across(all_of(group_cols_actual))) %>%
          dplyr::arrange(dplyr::desc(.data[[score_cols$fc]]), .data[[score_cols$padj]], .by_group = TRUE) %>%
          dplyr::slice_head(n = top_n_peaks) %>%
          dplyr::ungroup()
      } else {
        warning("No fold-change / adjusted-P columns detected in ", sheet, "; using all peaks.")
      }
      
      if (nrow(df2) > 0) {
        set_key <- do.call(paste, c(df2[group_cols_actual], sep = "__"))
        split_idx <- split(seq_len(nrow(df2)), set_key, drop = TRUE)
        
        sets <- lapply(names(split_idx), function(k) {
          idx <- split_idx[[k]]
          gr <- GRanges(
            seqnames = normalize_chr(df2[[coord_cols[1]]][idx]),
            ranges   = IRanges(
              start = as.integer(df2[[coord_cols[2]]][idx]),
              end   = as.integer(df2[[coord_cols[3]]][idx])
            )
          )
          GenomicRanges::reduce(gr)
        })
        names(sets) <- names(split_idx)
        sets <- sets[lengths(sets) > 0]
        if (length(sets) > 0) return(sets)
      }
    }
  }
  
  exclude <- intersect(names(df), c("seqnames", "start", "end", "chr", "chromosome", "peak_id", "id", "name"))
  candidate_cols <- setdiff(names(df), exclude)
  
  sets <- list()
  for (col in candidate_cols) {
    gr <- parse_peak_vector(df[[col]])
    if (length(gr) > 0) {
      sets[[col]] <- GenomicRanges::reduce(gr)
    }
  }
  
  sets
}

# ============================================================
# ---- build all_tasks: top-N peaks per analysis ----
# ============================================================

all_tasks <- list()

for (i in seq_along(sheet_rules)) {
  rule <- sheet_rules[[i]]
  xlsx_path <- file.path(input_dir, rule$workbook)
  
  cat(sprintf("[%d/%d] Reading %s / %s\n", i, length(sheet_rules), rule$workbook, rule$sheet))
  
  if (!file.exists(xlsx_path)) {
    warning("Workbook not found: ", xlsx_path)
    next
  }
  
  set_list <- tryCatch(
    extract_sets_from_sheet(
      xlsx_path,
      rule$sheet,
      group_cols = rule$group_cols,
      top_n_peaks = top_n_peaks
    ),
    error = function(e) {
      warning("Failed to parse sheet ", rule$sheet, ": ", conditionMessage(e))
      list()
    }
  )
  
  if (length(set_list) == 0) {
    warning("No sets found in ", rule$workbook, " / ", rule$sheet)
    next
  }
  
  for (set_name in names(set_list)) {
    all_tasks[[length(all_tasks) + 1]] <- list(
      workbook   = rule$workbook,
      sheet      = rule$sheet,
      set_name   = set_name,
      set_id     = clean_id(paste(rule$workbook, rule$sheet, set_name, sep = "__")),
      group_cols = paste(rule$group_cols, collapse = ","),
      set_gr     = set_list[[set_name]]
    )
  }
}

cat("Total sets to test: ", length(all_tasks), "\n", sep = "")
if (length(all_tasks) == 0) stop("No sets found to analyze.")

# ============================================================
# ---- HOMER inputs: top-N peaks per analysis, HOMER peak format ----
# ============================================================

homer_dir <- file.path(out_dir, "homer_inputs")
dir.create(homer_dir, recursive = TRUE, showWarnings = FALSE)

homer_genome      <- "hg38"    # must be installed locally (configureHomer.pl -install hg38)
homer_size        <- 200       # fixed 200bp window centered on each peak -- standard convention
# in the literature for TF motif finding (also HOMER's own default)
homer_threads     <- 4         # -p passed to findMotifsGenome.pl
homer_extra_args  <- "-mask"   # extra findMotifsGenome.pl flags

# HOMER peak file format (min 5 tab-separated columns, no header):
#   1. Unique Peak ID   2. chromosome   3. start   4. end   5. strand (+/- or 0/1)
write_homer_peak_file <- function(gr, set_id, out_path) {
  if (length(gr) == 0) {
    warning("Skipping HOMER peak file for '", set_id, "': no peaks after top-N selection")
    return(invisible(FALSE))
  }
  
  strand_chr <- as.character(GenomicRanges::strand(gr))
  strand_chr[strand_chr == "*"] <- "+"   # unstranded peaks default to "+"
  
  df <- data.frame(
    PeakID = paste0(set_id, "_peak_", seq_along(gr)),
    chr    = as.character(GenomicRanges::seqnames(gr)),
    start  = GenomicRanges::start(gr),
    end    = GenomicRanges::end(gr),
    strand = strand_chr,
    stringsAsFactors = FALSE
  )
  
  write.table(
    df, out_path,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  invisible(TRUE)
}

cat(
  "\nWriting HOMER peak files (top ", top_n_peaks, " peaks per analysis, ",
  length(all_tasks), " sets total) to: ", homer_dir, "\n", sep = ""
)

homer_cmds     <- character(0)
homer_manifest <- data.frame(
  set_id = character(0), n_peaks = integer(0), peak_file = character(0),
  stringsAsFactors = FALSE
)

for (task in all_tasks) {
  peak_file <- file.path(homer_dir, paste0(task$set_id, ".homerpeaks.txt"))
  ok <- write_homer_peak_file(task$set_gr, task$set_id, peak_file)
  
  if (isTRUE(ok)) {
    homer_out_subdir <- file.path(homer_dir, paste0(task$set_id, "_HOMER"))
    homer_cmds <- c(
      homer_cmds,
      sprintf(
        "findMotifsGenome.pl %s %s %s -size %s -p %d %s",
        shQuote(peak_file), homer_genome, shQuote(homer_out_subdir),
        homer_size, homer_threads, homer_extra_args
      )
    )
    homer_manifest <- rbind(
      homer_manifest,
      data.frame(
        set_id = task$set_id, n_peaks = length(task$set_gr), peak_file = peak_file,
        stringsAsFactors = FALSE
      )
    )
  }
}

writeLines(homer_cmds, file.path(homer_dir, "run_all_homer.sh"))
write.csv(homer_manifest, file.path(homer_dir, "homer_inputs_manifest.csv"), row.names = FALSE)

cat(
  "Wrote ", nrow(homer_manifest), " HOMER peak files + run_all_homer.sh",
  " (chmod +x it, or run: bash ", file.path(homer_dir, "run_all_homer.sh"), ")\n", sep = ""
)
cat(
  "Note: since sets are reduce()'d to merge overlapping regions, n_peaks per file",
  " may be slightly below ", top_n_peaks, " -- see homer_inputs_manifest.csv for exact counts.\n", sep = ""
)