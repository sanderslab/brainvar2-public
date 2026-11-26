# Bulk RNA-seq alignment resources

This directory contains scripts used to prepare a GRCh38 reference compatible
with the project's bulk and single-cell RNA-seq annotation release, create a
STAR genome index, build FASTQ manifests, and run STAR alignment. See
`run_STAR.sh` and `make_manifest_file.sh` for the command-line behavior and
required columns/paths.

## Reference and software notes

- The recorded alignment used STAR 2.7.10b_alpha_230301.
- The reference was based on GRCh38.p13 and GENCODE v41, with contigs selected
  to match the corresponding Broad/GTEx reference convention.
- `subset_broad_genome.py`, `contig_sleuthing.R`,
  `extract_gencode_contigs_to_subset_to.R`, and
  `subset_genome_to_GTEx_no_hla_decoy_contigs.py` support reference preparation
  and contig comparison.
- `make_manifest_file.sh` prepares a manifest from FASTQ names;
  `run_STAR.sh` runs alignment from a manifest.

The original work used project-controlled storage for FASTQs, references,
indexes, and aligned BAMs. Those data and locations are not distributed here.
Configure paths for an authorized execution environment and verify the
reference files and STAR index before running. This README records the
available scripts and historical reference settings; it does not establish
that a new run is identical to a prior production run.
