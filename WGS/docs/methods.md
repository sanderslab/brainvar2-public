# Methods

Paired-end whole-genome sequencing reads were processed with NVIDIA Parabricks
4.0.1-1. The `fq2bam` workflow performed BWA-MEM alignment to the
`grch38_full_plus_hs38d1_analysis_set` reference (NCBI assembly accession
`GCA_000001405.15`, with the hs38d1 decoy), mate fixing, Picard
2.18.2-compatible duplicate marking, and base-quality score recalibration
using dbSNP 138, Axiom Exome Plus, and Mills/1000 Genomes known-sites
resources. Optical duplicates used a 2,500-pixel distance. One whole-genome
germline gVCF per WGS library was produced with the Parabricks implementation
of GATK HaplotypeCaller using `StandardAnnotation` and `StandardHCAnnotation`
annotations and GQ bands from 10 through 90 in increments of 10.

The same whole-genome gVCF was supplied to each of 25 contig-specific
joint-calling records. Per-library gVCFs were imported into one GATK GenomicsDB
workspace per GRCh38 contig with a batch size of 50 and workspace
consolidation. GATK 4.2.6.1
GenotypeGVCFs performed cohort joint genotyping with StandardAnnotation and
AS_StandardAnnotation annotations. Per-contig raw calls were merged for VQSR
model training. Separate SNP and INDEL models used the resource roles, priors,
annotations, and tranches listed in `docs/reference.md`. INDEL and SNP models
were applied at truth-sensitivity thresholds of 99.0 and 99.8, respectively,
and the recalibrated per-contig VCFs were merged in reference contig order.

The joint-calling implementation is a focused derivative of nf-core/sarek
3.0.1. GATK 4.2.6.1 entered the production path at GenomicsDBImport; it did not
perform the original per-library HaplotypeCaller step. The production joint run
recorded Nextflow 23.10.1. See `docs/provenance.md` for workflow lineage and
production validation.
