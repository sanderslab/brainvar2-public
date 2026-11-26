# Description: 
# Perform variant quality control to remove low-quality and potential somatic variants. 
# Retain variants with allele frequency (AF) > 0.01 (1%) after all quality filters are applied. 
# Usage (Google Cloud Dataproc): # hailctl dataproc submit <your-dataproc-cluster> hail_comVarQC_single.py
# By Shan Dong 2024-05-06
# Note: check git repo https://github.com/stephansanders/quality_check for description

import hail as hl
hl.init()

def run_comVarQC(vcf_in_uri, vcf_out_uri, meta_uri, plink_out_uri, res_pca_sorce_uri):
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
    # keep "PASS" in filters field
    mt = mt.filter_rows(hl.len(mt.filters) == 0)
    # hard-filters suggested by GATK for germline short variant: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
    cond_row = ((mt.info.QD >= 2.0)
              & (mt.info.FS <= 60.0)
              & (mt.info.SOR <= 3.0)
              & (mt.info.MQ >= 40)
              & (mt.info.MQRankSum >= -12.5) 
              & (mt.info.ReadPosRankSum >= -8.0))
    mt = mt.filter_rows(cond_row, keep = True)
    # filter on depth (normal WGS coverage 30-40X)
    mt = mt.filter_entries( (mt.DP < 10) | (mt.DP > 200), keep = False) 
    # genotype QC
    ab = mt.AD[1] / hl.sum(mt.AD)
    filter_condition_ab = ((mt.GT.is_hom_ref() & (ab <= 0.1)) |
                            (mt.GT.is_het() & (ab >= 0.25) & (ab <= 0.75)) |
                            (mt.GT.is_hom_var() & (ab >= 0.9)))
    mt = mt.filter_entries(filter_condition_ab)
    # clean-up (remove AC = 0 loci) 
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    mt = mt.drop('variant_qc')
    # common variants wiht MAF cutoff of 1%
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AF[1] > 0.01, keep = True)
    # not so far from Hardy-Weinberg equilibrium as to suggest sequencing error
    mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 1e-6, keep = True)
    # clean-up (remove AC = 0 loci) 
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    mt = mt.drop('variant_qc')

    # write result
    hl.export_vcf(mt, vcf_out_uri)
    hl.export_plink(mt, plink_out_uri, ind_id = mt.s)
    # same commond as 'plink --vcf /path/to/file.vcf --make-bed --out sample --const-fid --keep-allele-order'

# inputs
meta_uri = 'gs://shanucsf/20240310_BV2_VCF_sample_list_keep.txt'
vcf_in_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.vcf.gz'
vcf_out_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.varQC.AF.vcf.gz'
plink_out_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.varQC.AF'
res_pca_sorce_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.varQC.AF.PCAscores.txt'

run_comVarQC(vcf_in_uri, vcf_out_uri, meta_uri, plink_out_uri,res_pca_sorce_uri)