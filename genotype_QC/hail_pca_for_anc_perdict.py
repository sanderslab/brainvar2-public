# Description: 
# prepare ancestry perdiction inputs using hail
# Usage (Google Cloud Dataproc): # hailctl dataproc submit <your-dataproc-cluster> hail_pca_for_anc_perdict.py
# By Shan Dong 2025-09-24
# Note: downstream ancestry prediction scripts and README at https://github.com/stephansanders/identity_check


import hail as hl
hl.init()

# inputs
build = 'GRCh38'
tag = "BrainVar2"
p5k_intervals = 'gs://hail-sandbox/wgs-dnms/IBD/purcell5k_grch38_liftover_2021-09-14.interval_list'
mt_uri = 'gs://brainvar2/result_check/brainvar2.022924.joint_called.recalibrated.varQC.gnomad.mt'
# result1 - merged filtered whole genome dataset
temp_p5k_mt_uri = 'gs://brainvar2/result_check/Ancestry/20250924_{}_p5k.mt'.format(tag)
# result2 - full samples list for IBD
sample_table_uri = 'gs://brainvar2/result_check/Ancestry/20250924_{}_relatedness_full_sample_list.txt'.format(tag)
# read matrixtable (after varQC)
mt = hl.read_matrix_table(mt_uri)
# load P5k interval list
interval_table = hl.import_locus_intervals(p5k_intervals, reference_genome=build)
# create filtered mt
mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]), keep = True)
# calls need to be unphased
mt = mt.select_entries(GT = hl.unphased_diploid_gt_index_call(mt.GT.n_alt_alleles()))
# repartition
mt = mt.repartition(50)
# Write/read temp dataset
mt.write(temp_p5k_mt_uri, overwrite=True)
# export list of samples involved
mt.cols().rename({'s' : 'Sample'}).export(sample_table_uri)



def run_pca (tag, p5k_mt_uri, temp_mt_uri, res_samples_uri, res_pca_sorce_uri,res_pca_eigenvalues_uri) :
    # load Purcell 5k-filtered dataset for pca and ancestry perdiction
    mt = hl.read_matrix_table(p5k_mt_uri)
    mt = mt.select_cols() 
    mt = mt.annotate_cols(Dataset = tag)
    # keep only entry fields that exist in autism dataset
    mt = mt.select_entries('GT')
    # load Purcell 5k-filtered version of the gnomad_hgdp_1kg_callset
    gnomad = hl.read_matrix_table('gs://hail-sandbox/PCA/gnomad_hgdp_1kg_P5k_GRCh38_2021-10-25.mt')
    gnomad = gnomad.annotate_cols(Dataset = 'gnomad_hgdp_1kg')
    # drop sample annotations that don't exist in autism dataset
    gnomad = gnomad.drop('bam_metrics', 'sample_qc', 'gnomad_sex_imputation', 'gnomad_population_inference', 
                        'gnomad_sample_qc_residuals', 'gnomad_sample_filters', 'gnomad_high_quality',
                        'gnomad_release', 'relatedness_inference', 'hgdp_tgp_meta', 'high_quality')
    # keep only entry fields that exist in autism dataset
    gnomad = gnomad.select_entries('GT')
    # now join them
    mt = mt.union_cols(gnomad)
    # write temp file
    mt = mt.drop('info')
    mt = mt.repartition(50)
    mt.write(temp_mt_uri, overwrite=True)
    # read merged file
    mt = hl.read_matrix_table(temp_mt_uri)
    # export sample list
    mt.cols().rename({'s' : 'Sample'}).export(res_samples_uri)
    # run PCA
    eigenvalues, scores, loadings = hl.hwe_normalized_pca(mt.GT, k=20, compute_loadings=True)
    scores.flatten().export(res_pca_sorce_uri)
    # export eigenvalues to a text file
    hl.Table.parallelize(
        hl.literal(list(enumerate(eigenvalues, start=1))).map(
            lambda x: hl.struct(PC=x[0], eigenvalue=x[1])
        )
    ).export(res_pca_eigenvalues_uri)

# inputs
tag = "BrainVar2"
temp_p5k_mt_uri = 'gs://brainvar2/result_check/Ancestry/20250924_{}_p5k.mt'.format(tag)
temp_mt_uri = 'gs://brainvar2/result_check/Ancestry/pca_merged_temp1.mt'
# Reuslts:
res_samples_uri = 'gs://brainvar2/result_check/Ancestry/20250924_brainVar2_pca_sample_list.txt'
res_pca_sorce_uri = 'gs://brainvar2/result_check/Ancestry/20250924_brainVar2_pca_scores.txt'
res_pca_eigenvalues_uri = 'gs://brainvar2/result_check/Ancestry/20250924_brainVar2_pca_eigenvalues.txt'

# run
run_pca(tag, temp_p5k_mt_uri, temp_mt_uri, res_samples_uri, res_pca_sorce_uri,res_pca_eigenvalues_uri)
