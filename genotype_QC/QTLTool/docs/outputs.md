# Outputs

| Output | Description |
| --- | --- |
| `*.bam.mbv.txt` | Native QTLtools candidate-level report for each RNA BAM. |
| `*.bam.mbv.log` | QTLtools log for each RNA BAM. |
| `mbv_all.tsv.gz` | Lossless concatenation of all native report rows. |
| `figures/per_alignment/*.svg` | Vector figures used for manual review and PowerPoint assembly. |
| `figures/per_alignment/*.png` | Raster copies convenient for PowerPoint. |
| `figures/perc_*_consistent.*` | Cohort concordance heatmaps. |
| `docs/decision_records/` | Approved rendered evidence and manual decision records for public release. |

Raw production outputs remain controlled. Only explicitly approved rendered
artifacts copied to `docs/decision_records/` are public. The collector output
deliberately has no `matched_donor`, `status`, or automated decision field.
