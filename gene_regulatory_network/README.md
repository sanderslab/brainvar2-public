# Gene regulatory network analysis

This folder contains two inference entry points, their shared helpers, figure scripts, and a variant motif disruption analysis:

- `grn/run_phase_grn.R`: infer phase-level TF-target relationships using peaks around a target-gene TSS.
- `coloc/run_locus_tf.R`: run a related model using candidate peaks around supplied loci.
- `visualization/`: figure scripts for SOX5 and AHI1 GRN results.
- `variant_motif_disruption/`: motifbreakR analysis for the AHI1 C>A variant.

Both pipelines perform motif matching, metacell RNA/ATAC retrieval, phase-level batch correction, overlapping group-lasso selection, post-selection inference, repeated sampling, Cauchy p-value combination, and coefficient-score reporting.

## Repository layout

```text
gene_regulatory_network/
  grn/
    run_phase_grn.R
    R/                         # data retrieval, motif matching, inference helpers
  coloc/
    run_locus_tf.R
    R/                         # GRN scoring helper
  visualization/
    GRN_network_plot.R
    SOX5_GRN.R
  variant_motif_disruption/
    motifbreak_analysis.R
  README.md
```

Run commands from the repository root. Results are written to `results/grn/` and `results/coloc/` by default.

## Requirements

- R packages: `ArchR`, `GenomicRanges`, `IRanges`, `Matrix`, `Seurat`, `irlba`, `biomaRt`, `BSgenome.Hsapiens.UCSC.hg38`, `Biostrings`, `JASPAR2024`, `TFBSTools`, `progress`, `grpregOverlap`
- Python 3 packages: `anndata`, `pandas`, `scipy`
- MEME Suite: `fimo` and `fasta-get-markov` available on `PATH`

## Inputs

The default layout is:

```text
data/
  atac/
    <ArchR project>.rds
    <ATAC metadata>.tsv
    *-ArrowFiles-<batch>.arrow
  JASPAR2024_CORE_pfms_meme/
    *.meme
explore/explore_results/metacells/
  metacells_period_<period>.h5ad
  mapping_p<period>_<batch>.csv
```

Required content:

- The ArchR project and Arrow files must contain a `PeakMatrix` using hg38 coordinates.
- ATAC metadata must include `Sample`, `Period`, `rna_majorclass_v1`, and a cell-ID column (`cellNames`, `CellID`, `cell_id`, `Unnamed: 0`, or `X`).
- Each metacell H5AD must store metacells in observations, gene symbols in variables, and expression in `X`.
- Mapping CSVs must include `SEACell`, cell-type labels (`broad_cell_type` or `majorclass_v1`), and either `archr_cell_id_std` or both `library` and `cell_id_std`.
- Motifs must be individual MEME-format files.

Input locations can be configured without editing the code:

```bash
export ATAC_PROJECT_FILE=/path/to/project.rds
export ATAC_METADATA_FILE=/path/to/metadata.tsv
export ATAC_DIR=/path/to/arrow_files
export METACELL_DIR=/path/to/metacells
export MOTIF_DIR=/path/to/meme_motifs
```

Optional settings are `GENE_COORDINATES_FILE`, `GRN_RESULTS_DIR`, `COLOC_RESULTS_DIR`, `FIMO_BIN`, and `FASTA_GET_MARKOV_BIN`.

`GENE_COORDINATES_FILE` must be a CSV with `hgnc_symbol`, `chromosome_name`, `start_position`, `end_position`, and `strand`. If omitted, gene coordinates are queried from Ensembl.

## Phase GRN

```bash
Rscript grn/run_phase_grn.R BCL11B prenatal 100000 1e-4 0.05 0.5 20 0.05 ExN 200 0.05 10
```

Arguments, in order:

```text
target_gene phase region_bp motif_hit_p min_open_fraction peak_score_quantile
n_runs fdr_cutoff cell_type max_peaks fimo_scan_p n_cores
```

Defaults are `100000`, `1e-4`, `0.05`, `0.5`, `100`, `0.05`, `ALL`, `200`, `0.05`, and `10`. `prenatal` includes periods 3–6 and `postnatal` includes periods 8–15, matching the Methods.

## Locus-centered analysis

Provide an hg38 locus CSV:

```csv
rsid,chr,pos
rs6446275,chr3,49429438
```

Run:

```bash
Rscript coloc/run_locus_tf.R GPX1 loci.csv prenatal ExN 50000 5000 0.01 1e-4 0.05 50 10 20 0.05
```

Arguments after `cell_type`, in order:

```text
window_bp enhancer_proximity_bp min_open_fraction motif_hit_p fimo_scan_p
max_peaks n_cores n_runs fdr_cutoff
```

Defaults are `50000`, `5000`, `0.01`, `1e-4`, `0.05`, `50`, `10`, `20`, and `0.05`.

## Main outputs

- `*_significant_tfs.csv`: significant TFs with mean score, median score, selection frequency, and BH-adjusted p-value.
- `*_variable_score_summary.csv`: statistics for every selected TF, CRE, and unpenalized period covariate.
- `*_bootstrap_variable_details.csv`: coefficient and post-selection p-value from each run.
- `*_combined_p_values.csv`: Cauchy-combined and BH-adjusted p-values.
- `*_X_matrix.csv`, `*_Y_matrix.csv`, `*_Z_matrix.csv`, and `*_motif_matrix.csv`: cached model inputs.

Completed GRN stages are reused automatically when their parameters match the requested analysis.

## Methods cross-check

The phase pipeline implements the described metacell RNA/ATAC aggregation, peak filtering by median open-cell fraction and phase accessibility score, motif-hit matrix construction, target-gene exclusion from TF predictors, RNA RPCA and ATAC LSI batch correction, unpenalized period indicators, overlapping group-lasso selection, post-selection OLS, repeated runs, Cauchy p-value combination, and BH adjustment.

Items to reconcile with the Methods text:

- The scripts use a `JASPAR2024_CORE_pfms_meme` default motif directory; Methods specifies JASPAR 2026. The `MOTIF_DIR` environment variable can point to another directory, so confirm which motif release produced the reported results.
- Both entry points assign prenatal periods 3–6, matching the Methods.
- Phase GRN defaults to a maximum of 200 peaks, matching the text. The locus-centered entry point defaults to a maximum of 50 peaks, and uses a 0.01 minimum open-cell fraction; document these as locus-analysis-specific settings if that is what was used.
- The locus-centered pipeline performs RNA RPCA and ATAC LSI correction and uses the same inference helper, but it is a separate supplied-locus analysis rather than the phase TSS-window analysis. Its defaults are a 50 kb window and 5 kb enhancer proximity.
