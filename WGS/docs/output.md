# Outputs

## Parabricks stage

`results/parabricks/cram/` contains, per sample:

- `<sample>.deduped.sorted.cram` and `.crai`;
- `<sample>.recal.table`;
- `<sample>.duplicate_metrics.txt`.

`results/parabricks/gvcf/` contains:

- `<sample>.grch38.gatk.haplotypecaller.g.vcf.gz`;
- the matching `.tbi` index.

## Joint stage

`results/joint/final_callset/` contains:

- `<cohort>.jointcalled.recalibrated.GRCh38.vcf.gz`;
- the matching `.tbi` index.

`results/joint/vqsr_models/` contains the SNP and INDEL recalibration tables,
indexes, tranche files, and model plotting scripts emitted by GATK.

`results/pipeline_info/software_versions.yml` records the tools invoked by the
selected stage. Nextflow also writes execution report, timeline, trace, and DAG
files there for non-test runs.

GenomicsDB workspaces, raw per-contig VCFs, and intermediate recalibrated VCFs
are retained in the Nextflow work directory but are not copied to the final
results tree by default.
