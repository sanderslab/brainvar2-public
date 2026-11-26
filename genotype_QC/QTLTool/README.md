# brainvar2_mbv

Portable code for the BrainVar/BrainVar2 QTLtools MBV sample-matching analysis.
The repository preserves the workflow that produced the project results:

1. Run QTLtools MBV separately for every RNA alignment against the same
   whole-genome cohort VCF.
2. Concatenate the native candidate-level reports into one table.
3. Use a Jupyter notebook to create per-alignment plots and cohort heatmaps.
4. Organize the plots in the controlled project PowerPoint and identify pairings
   manually.

The code does not automatically assign donors or replace the manual review.

## Environment

The recommended environment contains QTLtools 1.3.1 and the notebook packages:

```bash
conda env create -f environment.yml
conda activate brainvar2_mbv
```

The QTLtools 1.3.1 image historically named in project code remains available
as `naotokubota/qtltools:1.3.1`. Its verified manifest digest is
`sha256:e84ebf8ef271fb6ec999c6aa230c0ea8280bb6ebd14f79c0827d7b00c3be805a`.

## 1. Run MBV

Pass one set of indexed RNA BAMs to the Bash loop. Shell globs expand normally:

```bash
bash bin/run_mbv.sh \
  --vcf /controlled/brainvar2.whole_genome.vcf.gz \
  --outdir /controlled/mbv/bulk_rna \
  /controlled/bulk_rna_alignments/*.bam
```

Run the script again with a different output directory for another alignment
set, matching the way the BrainVar and BrainVar2 alignments were handled.

The default command is equivalent to:

```bash
QTLtools mbv \
  --bam RNA_ALIGNMENT.bam \
  --vcf WHOLE_GENOME_COHORT.vcf.gz \
  --out RNA_ALIGNMENT.bam.mbv.txt \
  --log RNA_ALIGNMENT.bam.mbv.log \
  --filter-mapping-quality 15
```

The mapping-quality setting matches the recovered BrainVar RNA workflow. All
other biological filters retain QTLtools 1.3.1 defaults. The script prints each
command before running it, and QTLtools writes one report and log per BAM directly
to the requested output directory.

## 2. Collect reports

Concatenate all native report rows without assigning or ranking donors:

```bash
python bin/collect_mbv_reports.py \
  --reports /controlled/mbv/bulk_rna/*.mbv.txt \
  --output /controlled/mbv/bulk_rna/mbv_all.tsv.gz
```

The collector requires the same candidate donor set in every report by default.
The resulting table contains `rna_alignment`, `candidate_donor`, and every
native QTLtools measurement.

## 3. Plot and review

Open [plot_mbv_results.ipynb](notebooks/plot_mbv_results.ipynb), set `DATA` to
the collected table, and run the notebook. The paired
[plot_mbv_results.py](notebooks/plot_mbv_results.py) is the maintainable
Jupytext source. The notebook writes:

- one SVG and PNG scatterplot per RNA alignment;
- heterozygous and homozygous concordance heatmaps.

The strongest candidates are labeled for visual review. The notebook does not
select a pairing. Figures and manual decisions remain in the approved
controlled project environment.

See [manual review](docs/manual_review.md) for the boundary between generated
figures and the PowerPoint-based pairing record.

## Historical decision record

The repository includes the approved master documentation used for the
BrainVar/BrainVar2 sample-matching review:

[Joseph_022124 Brainvar2 Sample Matching - documenting.pdf](docs/decision_records/Joseph_022124%20Brainvar2%20Sample%20Matching%20-%20documenting.pdf)

The 212-page deck preserves the analysis overview, non-synthetic MBV plots,
cohort-level conclusions, and individual-case reasoning used during manual
pairing review. See the [decision-record inventory](docs/decision_records/README.md)
for scope, dates, and checksum information. Additional approved rendered figures
may be placed under `docs/decision_records/figures/`; raw reports and collected
candidate tables remain excluded.

## Historical validation

The retained bulk archive contains 261 non-empty native reports, each with 271
candidate donors. The new collector reads all 261 reports into exactly 70,731
candidate rows. The retained historical primary summary contains 255 consensus
matches and 6 unresolved cases. The untouched legacy consolidator’s secondary
rank fields are known to be incorrect and are not reproduced here.

See [historical provenance](docs/historical_provenance.md) and
[verification](docs/verification.md) for the evidence boundary.

## Tests

```bash
bash scripts/check_release.sh
```

The tests use synthetic text reports only. No BrainVar BAM, VCF, native report,
or candidate-level table is release content. The explicitly approved historical
decision PDF and any curated figures under `docs/decision_records/` are the only
non-synthetic analysis artifacts intended for public release.

## Release status

This is a `0.1.0` release candidate. Repository ownership, full author order,
citation metadata, and DOI require team approval before `1.0.0`.
