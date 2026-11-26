# Identity and sex QC

Identity and sex checks were performed as companion analyses and are not part
of the small-variant Nextflow DAG.

## RNA-DNA identity matching

BrainVar RNA-DNA matching used QTLtools MBV. Retained evidence includes
hundreds of non-empty MBV reports and consolidated comparison tables,
demonstrating successful execution. The study methods identify QTLtools
v1.3.1. The retained run records do not provide a more specific executable
fingerprint, so this repository does not assert one.

The sanitized rerun wrapper preserves QTLtools' native report and does not
interpret or rank candidate donors:

```bash
bash bin/run_qtltools_mbv.sh \
  --bam PRIVATE_RNA_ALIGNMENT.bam \
  --vcf PRIVATE_WGS_GENOTYPES.vcf.gz \
  --out PRIVATE_SAMPLE.mbv.txt \
  --log PRIVATE_SAMPLE.mbv.log
```

Use `--include-samples`, `--exclude-samples`, and `--reg` when required by the
controlled analysis design. The BAM/SAM and VCF/BCF must be coordinate
compatible and indexed as required by QTLtools.

The legacy consolidation script is outside this software release. The public
wrapper preserves the native MBV report and does not apply a secondary ranking
or donor-classification rule. Any reported match classification must therefore
define its rule separately from this wrapper.

## Sex checks

Sex checking compared bulk-RNA XIST coverage with WGS chromosome-Y coverage.
Coverage was collected with `samtools coverage`; the
XIST interval used was `chrX:73824624-73830651` on GRCh38, and chromosome Y was
queried as `chrY`.

```bash
bash bin/collect_region_coverage.sh \
  --alignment PRIVATE_RNA_ALIGNMENT.bam \
  --region chrX:73824624-73830651 \
  --out PRIVATE_RNA_XIST.coverage.tsv

bash bin/collect_region_coverage.sh \
  --alignment PRIVATE_WGS_ALIGNMENT.cram \
  --reference Homo_sapiens_assembly38.fasta \
  --region chrY \
  --out PRIVATE_WGS_chrY.coverage.tsv
```

The study describes bimodal clustering of these measurements. This software
release provides the documented coverage extraction but does not implement a
classification threshold or clustering model. Any reported sex classification
must define that downstream rule separately.

## Evidence boundary

Native MBV and coverage outputs contain controlled sample-level information.
They belong in the approved data environment, not this public software
repository. Public evidence should be limited to deidentified aggregate
counts, methods, and figures approved by the study team.
