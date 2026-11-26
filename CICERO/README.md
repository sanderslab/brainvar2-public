# CICERO co-accessibility workflow

The numbered files at the folder root are the primary workflow corresponding
to the CICERO Methods. 

## Numbered Methods workflow

1. **Export grouped accessibility data:** run `Rscript 01_export_for_cicero.R`
   in the ArchR environment. It keeps groups with at least 500 cells and
   samples at most 5,000 cells per group (seed 123), then writes
   `cicero_export/`.
2. **Run CICERO:** make `cicero_export/` run `Rscript 02_run_cicero_from_export.R`. The script uses `k = 50`, a
   1 Mb window, a 500 kb distance constraint, 100 sampled windows, and an
   absolute co-accessibility threshold of 0.2.
3. **Shared rerun helper:** `03_cicero_core_pipeline.R` implements the same
   CICERO steps for group-specific reruns. It is sourced by step 4.
4. **Rerun the four specified groups:** run
   `Rscript 04_rerun_group_specific_k.R` with `--k_override` set to the
   Methods values: `Micro_prenatal:10,SMC_postnatal:25,Immune_postnatal:34,VLMC_postnatal:39`.
   The script's default ratio-based calculation does not produce all those
   exact values.
5. **Assemble the final result set:** run
   `./05_assemble_corrected_dataset.sh <original_dir> <corrected_dir>
   <merged_outdir>` to combine the rerun groups with the unchanged groups.

## Other code

- `exploratory/` contains QC summaries, cross-cell-type and developmental
  comparisons, sensitivity analyses, metadata reporting, and plotting helpers.
  These are supporting analyses, not steps in the Methods workflow above.
- `diagnostics/` contains checks of cell-count effects and the fixed `k`
  setting.

