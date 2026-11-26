# FGSEA trait enrichment

The retained analysis scripts are:

- `1.RNA_fgsea_all_markers.R` — RNA marker enrichment for disorder-associated gene sets.
- `2.ATACT_fgsea_gene_activity_v2.R` — ATAC gene-activity marker enrichment.

Both scripts include multiple pipeline branches. For the marker-file workflows used for cell-type-predominant and cell-type-by-development marker sets, the scripts first retain genes detected in at least 20% of cells in either comparison group, then rank those genes by signed `-log10(adjusted P)` and run `fgseaMultilevel` with gene-set sizes 5–500 and `eps = 0`. These workflows do not apply an additional adjusted-P significance cutoff before ranking.
