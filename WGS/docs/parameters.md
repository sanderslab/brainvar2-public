# Reference parameters

The following joint-calling parameters are required in addition to the common
parameters documented in `docs/usage.md`:

| Parameter | Required file |
| --- | --- |
| `--fasta` | GRCh38 reference FASTA |
| `--fai` | FASTA index |
| `--dict` | sequence dictionary |
| `--dbsnp` | dbSNP 138 VCF |
| `--dbsnp_index` | dbSNP tabix index |
| `--hapmap` | HapMap 3.3 VCF |
| `--hapmap_index` | HapMap tabix index |
| `--omni` | Omni 2.5 VCF |
| `--omni_index` | Omni tabix index |
| `--one_thousand_g` | 1000 Genomes phase 1 high-confidence SNP VCF |
| `--one_thousand_g_index` | 1000 Genomes tabix index |
| `--axiom` | Axiom Exome Plus polymorphism VCF |
| `--axiom_index` | Axiom tabix index |
| `--mills` | Mills/1000 Genomes gold-standard INDEL VCF |
| `--mills_index` | Mills tabix index |

Use `conf/resources.example.config` as the canonical fill-in template. All
files must conform to the same GRCh38 sequence dictionary.
