# GRCh38 reference contract

The workflow is intentionally restricted to the production
`grch38_full_plus_hs38d1_analysis_set` resource contract. The FASTA derives
from NCBI assembly accession `GCA_000001405.15` (GRCh38) with the hs38d1 decoy,
under the UCSC-ID alignment-pipeline distribution. The large reference files
are not stored in the software repository. They are installed from public
iGenomes/Broad resources with exact size and SHA-256 validation:

```bash
bash scripts/utility/download_resources.sh
bash scripts/utility/download_resources.sh --verify-only
```

The canonical paths, direct public source, file sizes, and checksums are listed
in `resources/README.md`. The ready-to-use Nextflow configuration is
`conf/resources.downloaded.config`.

## Parabricks stage

The production reference directory was extracted from
`grch38_full_plus_hs38d1_analysis_set_bwamem1.tar.gz` and contains
`Homo_sapiens_assembly38.fasta` plus its Parabricks-compatible BWA index files.
The public downloader reconstructs the equivalent layout by downloading the
FASTA and `.64.*` indexes from iGenomes, then linking the FASTA, FAI, and
dictionary into the BWA-index directory.
The BQSR known-sites inputs are:

- dbSNP build 138;
- Axiom Exome Plus population polymorphisms;
- Mills and 1000 Genomes gold-standard INDELs.

## Joint stage

The FASTA, `.fai`, sequence dictionary, VCFs, and tabix indexes must use the
same contig names and sequence dictionary. The production VQSR resource roles
are:

| Mode | Resource label | Known | Training | Truth | Prior |
| --- | --- | ---: | ---: | ---: | ---: |
| SNP | HapMap 3.3 | false | true | true | 15.0 |
| SNP | Omni 2.5 | false | true | true | 12.0 |
| SNP | 1000 Genomes phase 1 high-confidence SNPs | false | true | false | 10.0 |
| SNP | dbSNP 138 | true | false | false | 7.0 |
| INDEL | Axiom Exome Plus | false | true | false | 10.0 |
| INDEL | Mills and 1000 Genomes gold-standard INDELs | false | true | true | 12.0 |
| INDEL | dbSNP 138 | true | false | false | 2.0 |

The workflow stages every resource with its index and builds the GATK resource
labels from this fixed policy. A resource substitution is a methodological
change and should be documented and revalidated before release.

The public `dbsnp_138.hg38.vcf.gz` has different compressed bytes from the
local production copy, but their uncompressed VCF content has the same SHA-256:
`0ff368a3a4e16fc539be19ddbc58996d7e9a959a50e8597a7de07e217cce8527`.
The downloader uses the public gzip and its matching public tabix index.

Protected FASTQs are outside this resource contract. They remain user-supplied
controlled-access inputs and are never fetched or archived by the repository.

## Container identities

| Component | Immutable image reference |
| --- | --- |
| Parabricks | `nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e` |
| GATK | `quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db` |

The Parabricks image's anonymous remote availability and hardware contract are
documented in `docs/execution_environment.md`.
