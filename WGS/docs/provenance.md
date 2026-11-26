# Workflow provenance

## Production relationship

This repository is the maintained, portable implementation of the BrainVar2
GRCh38 germline small-variant workflow. The production callset predates this
standalone software repository. The release was assembled from the retained
production script, final joint-calling manifest, run configuration,
software-version records, and callsets; this Git history is therefore not
itself the historical execution record.

The released implementation preserves the production command path, tool
versions, reference contract, resource roles, annotations, thresholds, and
sample-to-contig orchestration. Direct comparisons with retained production
artifacts are recorded in [verification.md](verification.md).

The retained evidence establishes that:

- the final manifest contains 6,775 rows representing 271 WGS libraries across
  25 contigs;
- each manifest row for a library points to the same indexed whole-genome gVCF,
  rather than an independently called per-contig gVCF;
- all libraries used the same Parabricks `fq2bam` and HaplotypeCaller command
  path;
- production per-library jobs requested two GPUs, 24 CPU cores, 200 GB RAM,
  and 400 GB disk, with additional RAM or disk required for some libraries;
- two retained gVCF storage roots reflect storage and sample-name history, not
  different variant-calling methods; and
- the cohort stage used the nf-core/sarek-derived GenomicsDBImport,
  GenotypeGVCFs, and VQSR implementation encoded here.

The public repository contains no production manifests, sample identifiers,
machine-local paths, individual-level genomic files, or final cohort VCFs.
Cohort terminology and metadata versions are documented in
[cohort_accounting.md](cohort_accounting.md).

## Tool and workflow versions

| Component | Production version or revision |
| --- | --- |
| NVIDIA Parabricks image | `nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1` |
| NVIDIA Parabricks digest | `sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e` |
| Nextflow | 23.10.1 |
| nf-core/sarek base | 3.0.1 |
| GATK | 4.2.6.1 |
| GATK container digest | `sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db` |
| nf-core/sarek source tag and base commit | 3.0.1, `ad2b34f3` |
| Retained derivative source revision | `c01e5a2d` |

## Release scope

The source release is intentionally limited to the production GRCh38
small-variant workflow and its dependency closure. Experimental callers,
alternative references, exploratory workflows, caches, and controlled fixtures
are outside this repository's scientific scope.

The GATK modules preserve the production 4.2.6.1 tool version and operational
arguments while replacing machine-specific paths and implicit sample-count
parameters with validated channels. Large intermediates are not published by
default.

## Validation scope

One complete WGS library reproduced its retained indexed gVCF byte-for-byte.
For joint calling, all 271 retained gVCFs were rerun on GRCh38 chr21. The
regenerated and retained callsets had identical sample order, 654,028 variant
keys, genotype calls, and complete per-sample FORMAT payloads.

The retained production chr21 file carries filters from a genome-wide VQSR
model. The corresponding historical genome-wide recalibration tables and
tranche files were not retained as validation inputs. The validation therefore
establishes exact raw joint-genotype equivalence and successful execution of
the encoded VQSR path, but does not claim byte-identical regeneration of the
historical genome-wide VQSR model.
