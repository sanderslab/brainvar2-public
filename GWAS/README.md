# GWAS gene analysis with FUMA and MAGMA

This guide documents the FUMA SNP2GENE steps and MAGMA output used by `MAGMA_gwas_qq_plots.R`. The local R script reads existing MAGMA gene results and creates plots.


## FUMA website workflow

1. Prepare the full, harmonized GWAS summary-statistic file. Use all eligible SNPs for MAGMA gene analysis, not only genome-wide significant hits. FUMA requires a P-value column and either rsID or chromosome plus GRCh37 position; provide effect and non-effect allele columns when available. Confirm the input build and column names before upload.
2. Open FUMA’s **SNP2GENE** module and upload the summary-statistic file. Set the input genome build and map the P, chromosome, position, rsID, and allele columns to the matching fields in the file.
3. Select the population-matched reference panel used for that GWAS. The chosen reference panel affects LD correction in MAGMA; record the panel for each trait.
4. In the gene annotation settings, choose protein-coding genes and the gene annotation version used for the analysis. Set the MAGMA gene window to **1 kb upstream and 1 kb downstream**.
5. Enable MAGMA analysis. Use the standard MAGMA gene test settings reported for FUMA unless the original analysis records specify otherwise. The Methods do not report a gene-property or gene-set analysis as part of this step.
6. Submit the job, wait for completion, and download the MAGMA gene analysis results. Keep the job settings, FUMA/MAGMA versions, selected reference panel, and downloaded outputs with the analysis records.

FUMA’s tutorial describes SNP2GENE input requirements and MAGMA settings; the official [FUMA tutorial](https://fuma.ctglab.nl/tutorial) and [SNP2GENE interface](https://fuma.ctglab.nl/) are the references for the current website.

## Local MAGMA output and plots

`MAGMA_gwas_qq_plots.R` expects one `magma.genes.out` file per trait, with this layout:

`<GWAS_MAGMA_DIRECTORY>/<TRAIT>/magma.genes.out`

The script combines the gene-level P values, summarizes counts at several nominal thresholds, and creates gene-level P-value distribution and Q–Q plots. The Q–Q calculation excludes missing and zero P values, as described in the Methods. 

## FUMA input references

- [Official FUMA tutorial: SNP2GENE input and MAGMA analysis](https://fuma.ctglab.nl/tutorial)
- [Official FUMA SNP2GENE page](https://fuma.ctglab.nl/)
