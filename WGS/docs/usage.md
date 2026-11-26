# Usage

## Execution stages

`--stage parabricks` converts paired FASTQs to a recalibrated CRAM and gVCF for
each sample. `--stage joint` consumes indexed gVCFs and creates a recalibrated
cohort VCF. The stages are separate because the production run persisted and
transferred gVCFs between computing environments.

The production FASTQs and their manifest are protected inputs. Keep both in the
approved controlled-access environment; never commit them or add them to a
public source repository or software archive. The tracked FASTQs under
`tests/fixtures/fastq` are tiny synthetic topology fixtures.

## Public resource acquisition

Everything else needed by the workflow is publicly obtainable. From the
repository root:

```bash
bash scripts/utility/download_resources.sh
bash scripts/utility/download_resources.sh --verify-only
```

Run either stage with `-c conf/resources.downloaded.config` to use the
verified repo-local bundle. See `resources/README.md` for direct source URLs,
canonical paths, sizes, checksums, and the dbSNP compression-equivalence note.
Use an ignored site-specific `conf/resources.config` only when the same
resources are already managed elsewhere.

## Common parameters

| Parameter | Meaning |
| --- | --- |
| `--stage` | Required: `parabricks` or `joint` |
| `--input` | Required stage-specific CSV manifest |
| `--outdir` | Output directory; default `results` |
| `--cohort` | Cohort prefix for joint outputs; default `brainvar2` |
| `--contigs` | Ordered comma-separated contig list |

The production contig order is `chr1` through `chr22`, followed by `chrX`,
`chrY`, and `chrM`. If a subset is supplied, the joint samplesheet must contain
the complete sample-by-subset matrix.

## Parabricks parameters

| Parameter | Meaning |
| --- | --- |
| `--parabricks_refdir` | Directory containing the FASTA and Parabricks BWA index files |
| `--fasta_basename` | FASTA basename inside the reference directory |
| `--known_sites` | Comma-separated dbSNP, Axiom, and Mills known-sites VCFs |
| `--known_sites_indexes` | Comma-separated tabix indexes in the same order as `--known_sites` |

Each manifest row represents one paired FASTQ chunk. Repeat the same sample
identifier on multiple rows to align all chunks together; an exact duplicate
sample-FASTQ row is rejected. Each R1 header supplies the flowcell and lane used
to reconstruct that chunk's production read-group format. Sample and library
fields use the manifest sample identifier, and the platform is `illumina`.
The production manifest is protected because it records sample identifiers and
FASTQ locations.

## Joint-calling parameters

The joint stage requires `--fasta`, `--fai`, `--dict`, `--dbsnp`, and
`--dbsnp_index`, plus indexed HapMap, Omni, 1000 Genomes, Axiom, and Mills VCFs.
See `conf/resources.example.config` for every parameter and
`docs/reference.md` for the role and prior assigned to each resource.

## Profiles

- `standard`: local execution with tools already available on the host.
- `docker`: containerized GATK joint calling.
- `gpu_docker`: Docker with GPU passthrough for Parabricks.
- `singularity`: containerized execution with automatic mounts.
- `gpu_singularity`: Singularity with NVIDIA GPU passthrough for Parabricks.

Cluster executors and resource overrides belong in a site-specific config, not
in the source workflow. Each Parabricks process requests an NVIDIA GPU resource
by default; sites may override the GPU count and other resources according to
their scheduler and available hardware. The retained production request used
24 CPU threads, 200 GB host RAM, and 400 GB disk. Unusually large libraries may
require more RAM or disk. See `docs/execution_environment.md` for
version-matched NVIDIA requirements and the successful BrainVar2 execution
record.

`fq2bam` passes half of its task memory allocation to Parabricks via
`--memory-limit`. The remaining half is deliberate headroom for the wrapper,
sorting, CRAM compression, duplicate marking, and BQSR processes that share the
same task/container allocation.

## Restarting

Nextflow caches completed tasks. Add `-resume` to the original command after a
transient failure. Keep the matching `work/` directory until the run and its
output validation are complete.

## Companion QC and cohort utilities

The aggregate cohort reconciler and controlled-environment MBV/coverage
wrappers are documented in `docs/cohort_accounting.md` and
`docs/identity_qc.md`. They are not stages of the WGS workflow and do not alter
its biological outputs.
