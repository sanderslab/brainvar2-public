# Description: 
# Perform variant quality control to remove low-quality and potential somatic variants. 
# Generate trio table view
# Usage (Google Cloud Dataproc): # hailctl dataproc submit <your-dataproc-cluster> hail_get_VarQC_gnomadAnno_circos_input.py
# By Shan Dong 2025-08-19

import hail as hl
hl.init()

def run_varQC(vcf_in_uri, mt_out_uri, meta_uri):
    # read vcf
    mt = hl.import_vcf(vcf_in_uri, reference_genome='GRCh38', force_bgz=True, find_replace=('nul', '.'))
    mt = hl.split_multi_hts(mt)
    # filter low complexity regions
    lcr = hl.import_bed('gs://hail-sandbox/wgs-dnms/LCR-hs38.bed', reference_genome='GRCh38')
    mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)
    # filter centromeres
    ctm = hl.import_bed('gs://hail-sandbox/wgs-dnms/20231210.hg38.centromeres_UCSCTB.bed', reference_genome='GRCh38')
    mt = mt.filter_rows(hl.is_defined(ctm[mt.locus]), keep=False)
    # filter telomeres
    tlm = hl.import_bed('gs://hail-sandbox/wgs-dnms/20231201.hg38.telomeres_UCSCTB.bed', reference_genome='GRCh38')
    mt = mt.filter_rows(hl.is_defined(tlm[mt.locus]), keep=False)
    # samples to keep
    meta = hl.import_table(meta_uri).key_by('SampleID')
    mt = mt.annotate_cols(pheno=meta[mt.s])
    mt = mt.filter_cols(mt.pheno.Role != '', keep = True)
    # variant call rate filter
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95, keep = True)
    # only "PASS" in filter field
    mt = mt.filter_rows(hl.len(mt.filters) == 0)
    # hard-filters suggested by GATK for germline short variant: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    cond_row = ((mt.info.QD >= 2.0)
              & (mt.info.FS <= 60.0)
              & (mt.info.SOR <= 3.0)
              & (mt.info.MQ >= 40)
              & (mt.info.MQRankSum >= -12.5) 
              & (mt.info.ReadPosRankSum >= -8.0))
    mt = mt.filter_rows(cond_row, keep = True)
    # annotate with gnomAD genome frequencies (gnomad v4.1) https://discuss.gnomad.broadinstitute.org/t/difference-between-popmax-in-the-vcf-and-popmax-grpmaxfaf-in-the-website/63
    gnomad_ht4_genome = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht')
    mt = mt.annotate_rows(gnomadg_v4_max_AF = gnomad_ht4_genome.index(mt.row_key).grpmax.AF)
    # filter on depth (normal WGS coverage 30-40X)
    mt = mt.filter_entries( (mt.DP < 10) | (mt.DP > 200), keep = False) 
    # genotype QC
    ab = mt.AD[1] / hl.sum(mt.AD)
    filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                            (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                            (mt.GT.is_hom_var() & (ab >= 0.9)))
    mt = mt.filter_entries(filter_condition_ab)
    mt = mt.drop('variant_qc')
    mt = hl.variant_qc(mt)
    # not so far from Hardy-Weinberg equilibrium as to suggest sequencing error
    mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6, keep = True)
    # call rate > 0.95
    mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95, keep = True)
    # clean-up (remove AC = 0 loci) 
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    # write result mt
    mt.write(mt_out_uri, overwrite=True)

# inputs
meta_uri = 'gs://shanucsf/20240310_BV2_VCF_sample_list_keep.txt'
vcf_in_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.vcf.gz'
mt_out_uri = 'gs://brainvar2/result_check/brainvar2.022924.joint_called.recalibrated.varQC.gnomad.mt'
run_varQC(vcf_in_uri, mt_out_uri, meta_uri) #32036453 rows

# prepare ouput format for circos plot
# load hail matrixTable
mt = hl.read_matrix_table(mt_out_uri)
# extract row-level information
rows_table = mt.rows().key_by() 
rows_table.describe()
# select desired fields: chrom, position, AF (assuming AF is in the INFO field)
selected_table = rows_table.select(
    locus = rows_table.locus,
    alleles = rows_table.alleles,
    call_rate = rows_table.variant_qc.call_rate,
    variant_qc_AF = rows_table.variant_qc.AF[1],
    gnomadg_v4_max_AF = rows_table.gnomadg_v4_max_AF
)
# export to file
selected_table.export('gs://brainvar2/result_check/varINFO_AF_gnoamdAF.tsv.bgz')

