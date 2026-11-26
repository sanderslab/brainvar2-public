# LDSC heritability workflow

This folder contains the scripts and trait parameters used for the LDSC analyses. 
## Workflow

1. **Prepare trait-specific inputs** — `01_prepare_data/01_prepare_ldsc_data_by_trait.R` selects cohort SNPs shared with each GWAS summary-statistic dataset and the annotation data. It applies the documented MAF and summary-statistic QC filters, removes the MHC region, and writes trait-specific genotype, summary-statistic, and annotation inputs.
2. **Calculate LD scores** — `02_run_ldsc/02_calculate_ld_scores.R` generates baseline and annotation-specific LD scores for each trait.
3. **Estimate heritability** — `02_run_ldsc/03_run_heritability_ldsc.R` creates commands for the main heritability analysis, using trait prevalence parameters from `01_prepare_data/LDSC.parameters.xlsx` for binary traits.
4. **Run the overlap/coefficient analysis when needed** — `02_run_ldsc/04_run_heritability_ldsc_with_overlap.R` creates the alternate LDSC commands with annotation overlap and printed coefficients.
5. **Summarize outputs** — `03_summarize_results/05_summarize_heritability_results.R` parses the main LDSC logs and assembles summary tables. `03_summarize_results/06_summarize_overlap_results.R` parses overlap-analysis results and prepares selected summary tables.
6. **Plot results** — `04_plot_results/07_plot_ldsc_heatmap.R` creates the LDSC heatmap from the summarized results workbook.

