# Biological domain and NMF analysis

This folder currently contains the major-cell-type non-negative matrix factorization (NMF) analysis.

## Script

- `nmf_major_celltypes_only.py` builds a binary trait-by-gene matrix from significant RNA and ATAC FGSEA leading-edge results in Table S5, then decomposes it into five trait domains/programs and program gene-loading profiles.

The script selects major-cell-type specificity panels, retains significant trait-by-gene sets with at least 10 leading-edge genes, and runs 100 random NMF restarts with up to 1,000 iterations per run. It retains the solution with the lowest reconstruction error. Program stability is assessed by matching components across restarts with the Hungarian algorithm and cosine similarity. GO Biological Process enrichment is run on the top 25, 50, and 100 genes from each domain using the full NMF gene universe as background.

An MHC-excluded sensitivity analysis is run in addition to the MHC-included analysis. This allows comparison of the domains and gene universe with and without MHC-associated genes.

## Input and execution

The script requires an Excel workbook containing `Table S5` with trait, cell-type, stage, NES, adjusted-P, modality, and leading-edge gene fields. Run it with:

```bash
python nmf_major_celltypes_only.py \
  --workbook "/path/to/Supplementary Tables.xlsx" \
  --output ./NMF_major_celltypes_only
```

Python dependencies are listed in the script header: NumPy, pandas, SciPy, scikit-learn, Matplotlib, seaborn, requests, and openpyxl. GO enrichment uses the g:Profiler web API, so that step requires network access.

## Outputs

The output directory contains the binary trait-by-gene matrices, trait/program (`W`) and program/gene (`H`) loadings, consensus and restart-stability summaries, top-100 genes per program, GO enrichment tables, program labels, and a JSON file recording parameters and summary values. Results are written separately for MHC-included and MHC-excluded runs, with a comparison table at the output root.
