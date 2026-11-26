# Cell–cell communication analysis

This directory contains the scripts for donor-level CellChat inference,
differential communication testing between healthy and homozygous AHI1-risk
donors, and visualization of the results. Each donor is analyzed separately
to retain the donor as the unit of biological replication.

## Analysis workflow

The analysis uses the merged RNA H5AD object and its `lognorm` expression
layer. CellChat objects are generated for ExN, InN, Microglia, and OPC groups.
Within each donor, groups with fewer than 10 cells are omitted; donors with
fewer than two eligible groups are skipped.

Communication probabilities are estimated using a 10% truncated mean as the
primary model and triMean as a sensitivity model. Population-size weighting is
disabled and 100 resampling iterations are used. Differential analysis is
restricted to prenatal donors and compares healthy donors with homozygous
AHI1-risk donors. The model includes AHI1 status only; Sex and PMI are not
covariates.

## Numbered scripts

Run scripts from this directory, in order:

```bash
Rscript 01_run_cellchat_per_sample_h5ad.R
Rscript 02_differential_CCI_AHI1.R
Rscript 03_make_AHI1_plots.R
```

### Step 1: Per-donor CellChat

[01_run_cellchat_per_sample_h5ad.R](01_run_cellchat_per_sample_h5ad.R) reads
the H5AD `lognorm` layer, aligns the cell metadata to H5AD cell IDs, maps gene
IDs to symbols, and makes duplicate symbols unique. The expression matrix is
restricted to ligand and receptor genes from `CellChatDB.human`. The script
then creates per-donor CellChat objects, estimates ligand–receptor and pathway
communication, and saves run configuration and session information.

Outputs are written under `cellchat_results_01/`, including:

- `cellchat_objects/`: per-donor CellChat RDS objects, grouped by model and
  cell-group definition.
- `results/`: per-donor probability, p-value, pathway, and diagnostic tables.
- `gene_maps/`, `logs/`, and `config/`: gene mapping, run logs, sample lists,
  run configuration, and `sessionInfo.txt`.

The differential workflow uses the `truncatedMean10/subtype4` objects.
Broad-type and triMean objects are retained for sensitivity and visualization.

### Step 2: Differential communication

[02_differential_CCI_AHI1.R](02_differential_CCI_AHI1.R) selects prenatal
donors with CellChat objects and applies the following filters before testing:

- Ligand–receptor pairs must have nonzero communication in at least four
  donors.
- Directed ligand–receptor settings must have probability greater than 0.001
  in at least four donors.
- Each tested setting must have at least three non-missing donors in both the
  healthy and risk groups.
- Unavailable donor-setting measurements are retained as missing values.

Probabilities are bounded to [0, 1] and transformed with the arcsine
square-root transform. Linear models are fitted with `limma::lmFit()`,
moderated with `limma::eBayes(trend = TRUE, robust = TRUE)`, and adjusted for
multiple testing using the Benjamini–Hochberg procedure.

Results are written to
`cellchat_results_01/homo_differential/truncatedMean10/subtype4/prenatal/`.
This includes differential tables, detection and coverage summaries, the
sample metadata used, plot inputs, analysis configuration, and
`sessionInfo.txt`.

### Step 3: Figures and summary tables

[03_make_AHI1_plots.R](03_make_AHI1_plots.R) reads the step 1 CellChat objects
and step 2 results. It creates differential plots, pathway summaries,
per-sample visualizations, targeted glutamate/cholinergic summaries, and a
manifest of generated files under the step 2 results directory in `ALL_PLOTS/`.