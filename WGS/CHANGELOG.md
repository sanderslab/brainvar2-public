# Changelog

## 0.1.0-rc.1 - 2026-07-15

- Package the verified BrainVar2 GRCh38 production command path as a portable
  workflow.
- Add pinned Parabricks per-sample processing and GATK joint calling with VQSR.
- Add strict stage-specific manifest validation.
- Support multiple paired FASTQ chunks per sample with one read group per
  sequencing unit.
- Stage known-sites indexes explicitly and preserve container headroom by
  deriving Parabricks working memory from half of the task allocation.
- Add synthetic fixtures and complete stub topology checks.
- Document references, outputs, methods, provenance, and validation limits.
- Verify a complete three-chunk WGS sample and reproduce its retained pre-joint
  gVCF and index byte-for-byte.
- Verify 271-sample chr21 joint calling against production with identical
  sample order, 654,028 variant keys, genotype calls, and complete FORMAT
  payloads.
- Add a privacy-preserving joint-call comparison utility and publish the
  aggregate chr21 validation record.
- Verify anonymous access to the immutable Parabricks image and document its
  version-matched execution environment and BrainVar2-tested resource request.
- Require a GPU resource for Parabricks without fixing the portable workflow
  default to the two-GPU allocation used during validation.
- Document the BWA-MEM/GATK CPU-equivalent settings and indexed-gVCF handoff
  for users without Parabricks access.
