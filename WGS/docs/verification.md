# Validation record

Last updated: 2026-07-18

## Public automated validation

The consolidated repository check is:

```bash
bash scripts/check_release.sh
```

It validates:

- JSON schema and documented workflow parameters;
- positive and negative Parabricks and joint-manifest behavior;
- both complete workflow topologies with `nextflow run . -stub-run`;
- expected CRAM, gVCF, VQSR-model, and final cohort VCF output contracts;
- aggregate cohort-accounting logic with synthetic metadata;
- companion QTLtools MBV and regional-coverage command rendering;
- pinned container identities and production-critical arguments; and
- exclusion of controlled genomic files, machine-local paths, and credentials
  from the public source tree.

Continuous integration runs this check with Nextflow 23.10.1. Synthetic
fixtures validate software topology and contracts; they are not biological
benchmark data.

The pinned Parabricks 4.0.1-1 image remains available from NVIDIA NGC at digest
`sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e`.
An anonymous registry inspection reached its manifest, configuration, and
layers. A local GPU smoke test completed and reported `pbrun: 4.0.1-1`.

## Single-library pre-joint equivalence

One complete BrainVar2 WGS library with three paired FASTQ chunks was processed
through the released Parabricks path using the production GRCh38 reference and
dbSNP, Axiom, and Mills known-sites resources. The run used two NVIDIA L4 GPUs,
Nextflow 25.12.0-edge, and the pinned Parabricks 4.0.1-1 image.

- `fq2bam` completed in 1h02m05s with 121.9 GB peak RSS.
- HaplotypeCaller completed in 40m09s with 109 GB peak RSS.
- The CRAM passed `samtools quickcheck`, and its three header-derived read
  groups matched the retained production CRAM.
- The indexed gVCF contained one sample and 516,275,828 records.
- The regenerated gVCF and tabix index were byte-identical to their retained
  production counterparts.
- gVCF SHA-256:
  `4705cfa85d50f9adb4a053bb8bdbaeecc784f017c91826ed6acc1079de8bc526`.
- Index SHA-256:
  `b5a17802f6ec57cbb656675693e16887280dec58bbdf99c378dd7aa83c07e3ba`.

This comparison validates the complete per-library production-compatible path
through indexed gVCF output. Controlled reads and genomic outputs are not
distributed with the software repository.

## Complete-cohort chr21 joint-call equivalence

All 271 retained indexed gVCFs were processed through GenomicsDBImport and
GenotypeGVCFs for GRCh38 chr21 using Nextflow 23.10.1 and the pinned GATK
4.2.6.1 container. All eight workflow processes completed successfully.

The regenerated raw chr21 callset was compared with the retained production
chr21 file using [compare_joint_callsets.sh](../bin/compare_joint_callsets.sh).
Both callsets contained 271 samples and 654,028 records. They produced identical
aggregate values:

| Comparison | SHA-256 |
| --- | --- |
| Sample order | `01ff493eec2ff84664ae6caff8a08de4fc7ae3744292c7c29dc27768d1f3059e` |
| Variant keys plus every genotype call | `5534dea3e9f247c404dc15aec4ccfab40a2313882b3afd16e4f7336b00b1035e` |
| Variant keys plus complete per-sample FORMAT payloads | `a7986064b477df83345d3ea4e3da1453589859d8123a13aebd3c9c7a69dd5ecc` |

The machine-readable, sample-free result is
[chr21_joint_call_equivalence.tsv](../validation/chr21_joint_call_equivalence.tsv).
The comparison streams sample names and genotype rows directly through SHA-256;
it does not write or report individual identifiers or genotypes.

The FORMAT comparison excludes headers, QUAL, FILTER, and INFO. Those fields
are expected to differ because the retained production file carries the
original genome-wide VQSR annotations and filters. Record keys, genotype calls,
and every per-sample FORMAT value are identical.

## VQSR validation scope

The complete VQSR topology is covered by public stub tests, and the controlled
chr21 run successfully trained and applied new SNP and INDEL models using the
encoded GRCh38 resources, roles, annotations, priors, tranches, and
truth-sensitivity thresholds.

The retained production chr21 file was filtered by a model trained across the
whole genome. Its historical recalibration tables and tranche files were not
retained as comparison inputs. A newly trained chr21-only model is therefore a
functional test, not a byte-identical reproduction of the historical
genome-wide model. This limitation does not affect the exact raw joint-genotype
equivalence described above.
