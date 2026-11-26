# BrainVar2 WGS workflow

[![repository checks](https://github.com/JosephLalli/Brainvar2_wgs/actions/workflows/ci.yml/badge.svg)](https://github.com/JosephLalli/Brainvar2_wgs/actions/workflows/ci.yml)

This repository provides the portable, version-controlled implementation of
the BrainVar2 whole-genome germline small-variant workflow. It contains the
production GRCh38 path:

1. NVIDIA Parabricks 4.0.1-1 `fq2bam` for BWA-MEM alignment to the
   `grch38_full_plus_hs38d1_analysis_set` reference (NCBI assembly
   `GCA_000001405.15`, with the hs38d1 decoy), duplicate marking, and
   base-quality score recalibration, followed by Parabricks `haplotypecaller`
   in genomic VCF mode.
2. GATK 4.2.6.1 cohort joint genotyping and SNP/INDEL variant quality score
   recalibration derived from nf-core/sarek 3.0.1.

The public implementation was assembled from retained production scripts,
configuration, manifests, software-version records, and callsets. It was then
validated directly against retained production outputs. See
[workflow provenance](docs/provenance.md) and the
[validation record](docs/verification.md).

## Quick start: synthetic validation

Requirements: Nextflow 23.10.1 or later, Bash, and Python 3. The stub checks do
not require genomic software, reference downloads, containers, or a GPU.

```bash
git clone https://github.com/JosephLalli/Brainvar2_wgs.git
cd Brainvar2_wgs

python3 bin/validate_samplesheet.py \
  parabricks tests/fixtures/parabricks_samplesheet.csv
python3 bin/validate_samplesheet.py \
  joint tests/fixtures/joint_samplesheet.csv --contigs chr1

nextflow run . -stub-run -c conf/test_parabricks.config
nextflow run . -stub-run -c conf/test.config
```

These commands exercise both complete dataflow topologies with synthetic
fixtures. Stub files verify workflow wiring and output publication; they are
not biological test data and must not be used to evaluate variant calls.

## Controlled-data execution

Download and verify the public GRCh38 reference bundle once:

```bash
bash scripts/utility/download_resources.sh
bash scripts/utility/download_resources.sh --verify-only
```

This installs the approximately 12.5 GB canonical layout documented in
`resources/README.md`. Protected FASTQs are not downloaded; they must be
supplied through the approved controlled-access environment.

Parabricks stage:

```bash
nextflow run . \
  -profile gpu_docker \
  -c conf/resources.downloaded.config \
  --stage parabricks \
  --input samplesheets/parabricks.csv \
  --outdir results
```

Joint-calling stage:

```bash
nextflow run . \
  -profile docker \
  -c conf/resources.downloaded.config \
  --stage joint \
  --input samplesheets/joint.csv \
  --cohort brainvar2 \
  --outdir results
```

The GPU stage requires a suitable NVIDIA GPU resource, a compatible container
runtime, and access to `nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1`. The
exact image is currently downloadable from NVIDIA and pinned by digest. See
the [execution environment](docs/execution_environment.md) for the
version-matched driver/runtime minimums and the BrainVar2-tested allocation.
The joint stage uses `quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0` and does
not require a GPU.

## CPU-equivalent preprocessing and calling

Parabricks is not required to use or recreate the joint-called BrainVar2
cohort VCF. The stable handoff is one indexed, whole-genome germline gVCF per
WGS library. NVIDIA identifies BWA-MEM 0.7.15 and GATK 4.2.0.0 as the
[compatible CPU software](https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/compatiblecpusoftwareversions.html)
for Parabricks 4.0.1. A CPU implementation should reproduce the following
settings before supplying its gVCFs to this repository's unchanged GATK
4.2.6.1 joint stage:

| Step | Production-compatible CPU settings |
| --- | --- |
| Alignment | BWA-MEM 0.7.15 against the exact `grch38_full_plus_hs38d1_analysis_set` FASTA (NCBI assembly `GCA_000001405.15`, with the hs38d1 decoy); use `-Y` and NVIDIA's deterministic CPU-equivalence batch setting `-K 10000000` |
| Read groups | One read group per FASTQ pair: `ID=flowcell.sample.lane`, `PU=flowcell.lane`, `SM=sample`, `LB=sample`, `DS=Homo_sapiens_assembly38.fasta`, and `PL=illumina`; combine all pairs for the same library before duplicate marking |
| Alignment processing | Coordinate sort; add/fix mate information including MC and MQ tags; use Picard 2.18.2-compatible duplicate marking and an optical-duplicate pixel distance of 2,500 |
| BQSR | GATK BaseRecalibrator with dbSNP 138, Axiom Exome Plus population polymorphisms, and Mills/1000 Genomes gold-standard indels, followed by ApplyBQSR |
| Germline calling | GATK HaplotypeCaller 4.2.0.0 over the whole genome in GVCF mode, with `StandardAnnotation`, `StandardHCAnnotation`, and GQ bands 10, 20, 30, 40, 50, 60, 70, 80, and 90 |
| Unchanged cohort stage | Index each whole-genome gVCF and supply it to `--stage joint`; the manifest repeats that same gVCF for each of the 25 requested contigs |

The essential CPU HaplotypeCaller invocation is:

```bash
gatk HaplotypeCaller \
  -R Homo_sapiens_assembly38.fasta \
  -I sample.bqsr.cram \
  -O sample.grch38.gatk.haplotypecaller.g.vcf.gz \
  -ERC GVCF \
  -G StandardAnnotation \
  -G StandardHCAnnotation \
  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 \
  -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
  --native-pair-hmm-threads 16
```

The archived NVIDIA
[`fq2bam`](https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/man_fq2bam.html)
and
[`haplotypecaller`](https://docs.nvidia.com/clara/parabricks/4.0.1/documentation/tooldocs/man_haplotypecaller.html)
pages provide the corresponding CPU command structure. Do not substitute
BWA-MEM2, `samtools markdup`, different known-sites resources, a different
reference FASTA, or per-contig HaplotypeCaller runs and describe the result as
the production-equivalent path. Diploid ploidy and the other GATK defaults were
retained; no project-specific mapping-quality threshold was added. NVIDIA
documents CPU compatibility, but this repository's byte-for-byte validation
used Parabricks; an independent CPU rerun has not been compared byte-for-byte.

If the same resources already exist at site-managed paths, copy
`conf/resources.example.config` to the ignored
`conf/resources.config` and replace its placeholders instead of downloading
another bundle.

## Software and protected-data boundary

The public software repository and its archival releases contain the workflow,
portable configs, documentation, validation utilities, checksums, and tiny
synthetic fixtures. They do not contain BrainVar FASTQs, production samplesheets,
BAMs/CRAMs, gVCFs, cohort VCFs, or the downloaded GRCh38 bundle.

FASTQs are protected individual-level human genomic data. Keep them and the
private manifest that names them in the approved controlled-access repository.
The public resource downloader acquires only the reference FASTA, BWA indexes,
and GATK resource VCFs required to run the software.

## Validation

Controlled-data validation covered both workflow boundaries:

- one complete three-chunk WGS library reproduced its retained indexed gVCF
  byte-for-byte; and
- all 271 retained indexed gVCFs reproduced the chr21 production sample order,
  654,028 variant keys, genotype calls, and complete per-sample FORMAT payloads.

The second comparison evaluates the raw joint-genotype content before VQSR-only
annotations and filters. Public synthetic tests validate topology and contracts;
they do not distribute controlled genomic data. Exact comparison definitions,
aggregate hashes, and limitations are documented in
[docs/verification.md](docs/verification.md).

## Data and code availability

The workflow source is available under the MIT License at
`https://github.com/JosephLalli/Brainvar2_wgs`. Manuscripts using this workflow
should cite the exact release tag or commit that was reviewed.

Individual-level BrainVar sequencing files, alignments, genotypes, and
production manifests are controlled human genomic data and are not distributed
through this software repository. Their access conditions are governed by the
study approvals and the manuscript's Data and Code Availability statement.
All public reference resources required by the workflow are independently
downloadable with recorded sizes and SHA-256 checksums.

## Input manifests

The Parabricks manifest has one paired-end sequencing chunk per row. Repeat a
sample identifier on additional rows when that sample has multiple FASTQ pairs;
the workflow combines those chunks in one alignment task while assigning a
header-derived read group to each pair:

```csv
sample,fastq_1,fastq_2
sample_1,/data/sample_1_chunk1_R1.fastq.gz,/data/sample_1_chunk1_R2.fastq.gz
sample_1,/data/sample_1_chunk2_R1.fastq.gz,/data/sample_1_chunk2_R2.fastq.gz
```

This manifest contains protected file locations and sample identifiers. Do not
commit the production copy or include it in a software release.

The joint manifest is long-form, with one row per sample and requested contig.
The same whole-genome gVCF may therefore appear on multiple rows:

```csv
sample,contig,gvcf,gvcf_index
sample_1,chr1,/data/sample_1.g.vcf.gz,/data/sample_1.g.vcf.gz.tbi
sample_1,chr2,/data/sample_1.g.vcf.gz,/data/sample_1.g.vcf.gz.tbi
```

The workflow rejects duplicate sample-FASTQ records, duplicate joint records,
unexpected contigs, incomplete sample-contig matrices, missing paths, and
unsafe sample identifiers before compute-intensive tasks begin.

## Documentation

- [Usage and parameter details](docs/usage.md)
- [Reference parameter list](docs/parameters.md)
- [GRCh38 reference bundle](docs/reference.md)
- [Public resource acquisition and checksums](resources/README.md)
- [Outputs](docs/output.md)
- [Methods](docs/methods.md)
- [Production provenance](docs/provenance.md)
- [Cohort accounting](docs/cohort_accounting.md)
- [Identity and sex QC](docs/identity_qc.md)
- [Verification record](docs/verification.md)
- [Execution environment](docs/execution_environment.md)
- [Changelog](CHANGELOG.md)

## License

This derivative is distributed under the MIT License. See `LICENSE` and
`THIRD_PARTY_NOTICES.md` for the retained upstream attribution. The MIT terms
permit use, modification, and redistribution provided that the copyright and
permission notice is retained in copies or substantial portions of the code.
