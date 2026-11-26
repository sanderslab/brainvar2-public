# GRCh38 runtime resources

The BrainVar2 FASTQs are protected human genomic data and are never downloaded
or distributed by this repository. Everything else needed by the production
workflow is either tracked as source code and tiny synthetic test data or
installed from public resources.

## Automated acquisition

From the repository root, download and verify the approximately 12.5 GB
production-compatible GRCh38 bundle:

```bash
bash scripts/utility/download_resources.sh
bash scripts/utility/download_resources.sh --verify-only
```

The downloader requires Bash, `curl`, and `sha256sum`. It is idempotent:
verified files are skipped, interrupted `.part` downloads are resumed, and an
existing file that fails validation is never overwritten automatically.

The public source is the credential-free Illumina iGenomes AWS distribution of
the Broad GATK GRCh38 bundle:

```text
https://ngi-igenomes.s3.amazonaws.com/igenomes/Homo_sapiens/GATK/GRCh38/
```

The GATK resource-bundle documentation describes the reference and VQSR
resources:
<https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle>

## Canonical paths

```text
resources/GRCh38/
├── Sequence/
│   ├── WholeGenomeFasta/
│   │   ├── Homo_sapiens_assembly38.fasta
│   │   ├── Homo_sapiens_assembly38.fasta.fai
│   │   └── Homo_sapiens_assembly38.dict
│   └── BWAIndex/
│       ├── Homo_sapiens_assembly38.fasta -> ../WholeGenomeFasta/...
│       └── Homo_sapiens_assembly38.fasta.64.*
└── Annotation/GATKBundle/
    ├── dbsnp_138.hg38.vcf.gz
    ├── Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
    ├── Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
    ├── hapmap_3.3.hg38.vcf.gz
    ├── 1000G_omni2.5.hg38.vcf.gz
    └── 1000G_phase1.snps.high_confidence.hg38.vcf.gz
```

Each VCF has its matching `.tbi`. The downloader also links the FASTA, FAI,
and dictionary into `Sequence/BWAIndex` so Parabricks sees the reference and
its `.64.*` BWA indexes in one directory. Use
`conf/resources.downloaded.config` to select this layout.

## Verified manifest

Sizes are bytes. SHA-256 values identify the exact public objects downloaded by
the helper.

| Relative path under `resources/GRCh38` | Bytes | SHA-256 |
| --- | ---: | --- |
| `Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta` | 3249912778 | `93157a161863464c9435062fd67c173fdaf99cb8b32f1455018361387ffa5564` |
| `Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai` | 160928 | `edefd93c489dc1baefad312f40388089f8db5cf6dcc3ba0955669ead274e8b6b` |
| `Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict` | 581712 | `07523b4a0afd1127ae909f92f8b3911aff096461aed7eff3f1acf5bdfc1dc103` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.alt` | 487553 | `d9254da07b8030e26129dc29d9d02b9c30a360b233a367ce041691c00407d510` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.amb` | 20199 | `4199db9e913e372b2dae4ac17723def657f696daf71024b66f6b695947918c77` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.ann` | 455474 | `49231cd8462826afe424aeb20b6b2b9ecff32ff1a90f88f5d5f3c1a0cfb9ba99` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.bwt` | 3217347004 | `a01ecc6bfae7abcf3542f68314da843d5594f8b8cf3eb7799371ac3e99231754` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.pac` | 804336731 | `5dd6629f95e9765a19821fd1ffc3e883e5af177527a220a6b72a4e5e492114dc` |
| `Sequence/BWAIndex/Homo_sapiens_assembly38.fasta.64.sa` | 1608673512 | `b87d29d937b58dfde9943d8300ab0f9490bae5ba4ddc54726da023aeebd48da4` |
| `Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz` | 1560889937 | `04d74d585c2adbcc601e159d6c610d479598779749d7e67a26e35af6cf58a97c` |
| `Annotation/GATKBundle/dbsnp_138.hg38.vcf.gz.tbi` | 2321422 | `12c86c71249bd7a29059b213807e3681aadb3237c0ef60188fc144eabad7514f` |
| `Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz` | 3053999 | `8bca5d7af481f8bb0cf657d0665a6a84be7f0b24c254bf6a912c9733e39f4217` |
| `Annotation/GATKBundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi` | 421689 | `593660b27c1edca1201f2b779dec189837dee9dab993baf2e7b3200f87fc038c` |
| `Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz` | 20685880 | `e11195ef119927518f01965bfe63591ec4c761b727ed26e72f17fed23a6739e5` |
| `Annotation/GATKBundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi` | 1500013 | `6856f5ad5690e1f095960ecee450e26a67f9486a8336aa08be4f599d958c1bac` |
| `Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz` | 62043448 | `6c0402a3b77e97e1067f52fe1b13bbaf39eb4d8b846d24e76b3ea8fe82237188` |
| `Annotation/GATKBundle/hapmap_3.3.hg38.vcf.gz.tbi` | 1552123 | `2aed605c4b96aebd480d9ef515a59777387e17b6988423688b1993aadf2ea155` |
| `Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz` | 53238342 | `9c4aa8e66972591500cc46c2b6792a74d81ae856cc467e22f5f58169dc987399` |
| `Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi` | 1544108 | `60aaea8e7b62831671c4c129dbb1deb364e991afe185bf3a0a220916fb7e9cc3` |
| `Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz` | 1888262073 | `734d6634a98094ada4eb365a3a6aa4babcbecf3617a95fb5a23ed3ff98a595d1` |
| `Annotation/GATKBundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi` | 2128536 | `7abef1e31a3d513f393c0bccbfb80df12ab7eb141d175707e3dee87981a2f4bc` |

The production dbSNP gzip stream was not byte-identical to this public object.
Their uncompressed VCF content is identical, with SHA-256
`0ff368a3a4e16fc539be19ddbc58996d7e9a959a50e8597a7de07e217cce8527`.
The downloaded public gzip and its matching public tabix index are therefore
the portable pair used by `conf/resources.downloaded.config`.

## Data boundary

The downloaded reference tree is ignored by Git and is not part of the source
repository or an archived software release. Protected FASTQs and private production
samplesheets must be supplied through the approved controlled-access storage
environment. Do not copy either into this repository.
