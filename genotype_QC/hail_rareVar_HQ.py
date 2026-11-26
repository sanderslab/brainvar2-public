# Description: 
# Select high quality rare variants from WGS dataset
# Usage (Google Cloud Dataproc): # hailctl dataproc submit <your-dataproc-cluster> hail_rareVar_HQ.py
# By Shan Dong 2025-05-12
# Note: check git repo https://github.com/stephansanders/quality_check for description

import hail as hl
hl.init()

def get_rare_HQ_var(vcf_in_uri, vcf_out_uri, meta_uri, out_table_uri):
    # read vcf
    mt = hl.import_vcf(vcf_in_uri, reference_genome='GRCh38', force_bgz=True, find_replace=('nul', '.'))
    mt = hl.split_multi_hts(mt) # hail built-in split multi method
    # filter low complexity regions
    lcr = hl.import_bed('gs://hail-sandbox/wgs-dnms/LCR-hs38.bed', reference_genome='GRCh38')
    mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)
    # annotate with gnomAD genome frequencies (gnomad v4.1) https://discuss.gnomad.broadinstitute.org/t/difference-between-popmax-in-the-vcf-and-popmax-grpmaxfaf-in-the-website/63
    gnomad_ht4_genome = hl.read_table('gs://gcp-public-data--gnomad/release/4.1/ht/genomes/gnomad.genomes.v4.1.sites.ht')
    mt = mt.annotate_rows(gnomadg_v4_max_AF = gnomad_ht4_genome.index(mt.row_key).grpmax.AF)
    mt = mt.filter_rows( (mt.gnomadg_v4_max_AF <= 0.0001) | (hl.is_missing(mt.gnomadg_v4_max_AF)), keep = True)
    # samples to keep
    meta = hl.import_table(meta_uri).key_by('SampleID')
    mt = mt.annotate_cols(pheno=meta[mt.s])
    mt = mt.filter_cols(mt.pheno.Role != '', keep = True)
    # HQ filters (used for denovos -- loose filters)
    # filter for VQSR - PASS for SNVs only
    mt = mt.filter_rows((hl.is_snp(mt.alleles[0], mt.alleles[1]) &(mt.filters.size() == 0))
                    | hl.is_indel(mt.alleles[0], mt.alleles[1]))
    # filter on depth
    mt = mt.filter_entries( (mt.DP < 10) | (mt.DP > 200), keep = False) 
    # row (varinat INFO) level filters 
    dn_snv_cond_row = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                        & (mt.qual >= 200)
                        & (mt.info.SOR <= 2.5)
                        & (mt.info.ReadPosRankSum >= -1.4)
                        & (mt.info.QD >= 3.0)
                        & (mt.info.MQ >= 60))
    dn_indel_cond_row = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                        & (mt.qual >= 200)
                        & (mt.info.SOR <= 3)
                        & (mt.info.ReadPosRankSum >= -1.7)
                        & (mt.info.QD >= 4.0)
                        & (mt.info.MQ >= 50))
    mt = mt.filter_rows(dn_snv_cond_row | dn_indel_cond_row, keep = True)
    # homozygous filters
    ab = mt.AD[1]/hl.sum(mt.AD)
    hom_snv_cond = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_hom_ref()
                    & (mt.GQ >= 30.0)
                    & (ab <= 0.05))
    hom_indel_cond = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_hom_ref()
                    & (mt.DP >= 16)
                    & (mt.GQ >= 30.0)
                    & (ab <= 0.05))
    # heterzygous filters
    het_snv_cond = (hl.is_snp(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_het()
                    & (mt.GQ >= 99.0)
                    & (ab >= 0.22)
                    & (ab <= 0.78))
    het_indel_cond = (hl.is_indel(mt.alleles[0], mt.alleles[1])
                    & mt.GT.is_het()
                    & (mt.GQ >= 99.0)
                    & (ab >= 0.20)
                    & (ab <= 0.80))
    filter_condition = (hom_snv_cond | hom_indel_cond | het_snv_cond | het_indel_cond)
    mt = mt.filter_entries(filter_condition, keep = True)
    # GQ mean filter
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= 50, keep = True)
    # clean-up: remove AC = 0 loci
    mt = hl.variant_qc(mt)
    mt = mt.filter_rows(mt.variant_qc.AC[1] > 0, keep = True)
    mt = mt.drop('variant_qc')
    # export vcf
    hl.export_vcf(mt, vcf_out_uri)

vcf_in_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.vcf.gz'
meta_uri = 'gs://shanucsf/20240310_BV2_VCF_sample_list_keep.txt'
vcf_out_uri = 'gs://brainvar2/brainvar2.022924.joint_called.recalibrated.HQ.rare.vcf.gz'

# Run PCA 2025-07-01
res_pca_sorce_uri = 'gs://brainvar2/result_check/brainvar2.022924.joint_called.recalibrated.HQ.rare.PCAscores.txt'
# read vcf
mt = hl.import_vcf(vcf_in_uri, reference_genome='GRCh38', force_bgz=True, find_replace=('nul', '.'))
eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
scores.flatten().export(res_pca_sorce_uri)