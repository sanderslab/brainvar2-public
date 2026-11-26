# Methods

For each coordinate-sorted bulk RNA-seq BAM, QTLtools MBV was run against the
BrainVar whole-genome cohort genotype VCF. The RNA alignment and genotype VCF
used the same reference coordinate system. The recovered BrainVar RNA workflow
specified QTLtools 1.3.1 and minimum read mapping quality 15. Native QTLtools
defaults were used for the remaining MBV filters.

QTLtools produced one report per RNA alignment, with one row per candidate WGS
sample and measurements including covered heterozygous and homozygous genotypes
and their concordance with observed RNA alleles. The native reports were
concatenated into a long table without changing the QTLtools values.

Per-alignment scatterplots compared `perc_het_consistent` with
`perc_hom_consistent`; the strongest candidates were labeled. Cohort heatmaps
displayed each concordance metric across RNA alignments and candidate WGS
samples. Plots were manually organized in a controlled PowerPoint presentation,
and RNA-DNA pairings were identified by manual visual review. No automated
threshold or software classifier generated the final pairing decisions.

QTLtools should be cited as Delaneau et al., *Nature Communications* (2017),
DOI 10.1038/ncomms15452, PMID 28516912.
