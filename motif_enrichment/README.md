# HOMER motif enrichment

This folder contains the HOMER workflow for DAR and caQTL peak sets. Files are grouped by analysis stage.

## Workflow

### 1. Prepare peak inputs

- `01_prepare_inputs/01_prepare_dar_homer_inputs.R` reads the DAR result tables, selects peaks for each analysis group, writes HOMER peak files, and generates a DAR command list.
- `01_prepare_inputs/02_prepare_caqtl_homer_inputs.R` converts significant caQTL results into unique HOMER-format peak files by cell type and developmental stage.

### 2. Run HOMER

Run the generated DAR commands from `02_run_homer/run_homer_dars.sh` and the caQTL commands in `02_run_homer/run_homer_caqtl.sh` in an environment with HOMER and the hg38 genome installed. The command files contain the paths used when they were generated; update those paths if the inputs or output directories have moved.

### 3. Summarize and plot

- `03_summarize_plot/01_summarize_homer_results.R` reads HOMER result folders and produces integrated known and de novo motif tables, recurrence summaries, and overview outputs. Its `strong_enrichment` label is a supplementary summary based on q-value and a derived target/background percentage ratio; it does not change the HOMER results.
- `03_summarize_plot/02_plot_homer_motif_summaries.R` reads the integrated known-motif table, retains known motifs with Benjamini-adjusted *P* < 0.05 and log₂ odds ratio > 1, and creates motif significance/effect-size violin plots and summary visualizations. The log₂ odds ratio is calculated from HOMER target and background motif percentages.
