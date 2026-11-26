process Transpose {
    label 'Transpose'
    publishDir "${params.work_path}/${target_folder}", mode: 'copy', storageClass: 'INTELLIGENT_TIERING'
    // , tags: [TaskID: "${task.name.toString()}", Project: 'myproject']

    input:
    tuple val(sample_group), val(pc_name), val(cell_group), path(bed_file), val(sorted_bed_file), val(target_folder)
   
    output:
    tuple val(sample_group), val(pc_name), val(cell_group), emit: all_of_it
    path "${sorted_bed_file}", emit: sorted_bed_file

    script:
    """
    transpose.py ${bed_file} ${sorted_bed_file}
    """

    stub:
    """
    touch "${sorted_bed_file}"
    """
}

process SortBedFile {
    label 'SortBedFile'
    publishDir "${params.work_path}/${target_folder}", mode: 'copy'

    input:
    tuple val(sample_group), val(pc_name), val(cell_group), val(genotype), path(bed_file), val(sorted_bed_file), val(target_folder)
   
    output:
    tuple val(sample_group), val(pc_name), val(cell_group), val(genotype), emit: all_of_it
    path "${sorted_bed_file}", emit: sorted_bed_file

    script:
    """
    sort_bed_file.py ${bed_file} ${sorted_bed_file}
    """

    stub:
    """
    touch "${sorted_bed_file}"
    """
}

process FillNaNWithRowMedian {
    label 'FillNaNWithRowMedian'
    publishDir "${params.work_path}/${target_folder}", mode: 'copy'

    input:
    tuple val(sample_group), val(pc_name), val(cell_group), val(genotype), path(covariates), val(target_folder), val(target_file)

    output:
    tuple val(sample_group), val(pc_name), val(cell_group), val(genotype), emit: all_of_it
    path "${target_file}", emit: covariates_filled

    script:
    """
    fill_nan.py ${covariates} ${target_file}
    """

    stub:
    """
    touch ${target_file}
    """
}

process RunTensorqtl {
    label 'tensorqtl'
    publishDir "${params.work_path}/${target_folder}", mode: 'copy'

    input:
    tuple val(target_folder), val(file_name), val(genotype), path(plink_bed), path(plink_bim), path(plink_fam), path(covariates), path(phenotypes)
   
    output:
    path "${file_name}.tensorQTL.cis.log", emit: tensorqtl_log
    path "${file_name}.cis_qtl.txt.gz", emit: tensorqtl_file

    script:
    """
    python3 -m tensorqtl ${genotype} ${phenotypes} ${file_name} \
    --covariates ${covariates} \
    --permutations 1000 \
    --mode cis
    """

    stub:
    """
    touch "${file_name}.tensorQTL.cis.log"
    touch "${file_name}.cis_qtl.txt.gz"
    """
}

process RunTensorqtlNominal {
    label 'tensorqtl_nominal'
    publishDir "${params.work_path}/${target_folder}/cis_nominal", mode: 'copy'

    input:
    tuple val(target_folder), val(file_name), val(genotype), path(plink_bed), path(plink_bim), path(plink_fam), path(covariates), path(phenotypes)

    output:
    path "${file_name}.tensorQTL.cis_nominal.log", emit: tensorqtl_log
    path "${file_name}.cis_qtl_pairs.*.parquet", emit: tensorqtl_files

    script:
    """
    python3 -m tensorqtl ${genotype} ${phenotypes} ${file_name} \
    --covariates ${covariates} \
    --mode cis_nominal
    """

    stub:
    """
    touch "${file_name}.tensorQTL.cis_nominal.log"
    touch "${file_name}.cis_qtl_pairs.chr1.parquet"
    """
}

