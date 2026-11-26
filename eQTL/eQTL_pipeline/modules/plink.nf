process RunPlink {
    label 'plink'
    publishDir "${params.work_path}/plink_vcf", mode: 'copy', storageClass: 'INTELLIGENT_TIERING'

    input:
    tuple val(genotype), path(genotype_vcf), path(sample_information_file), val(plink_par)
   
    output:
    path "${genotype}.log", emit: plink_log
    path "${genotype}.bed", emit: genotype_bed
    path "${genotype}.bim", emit: variant_information
    path "${genotype}.fam", emit: sample_information

    script:
    """
    plink2 --vcf ${genotype_vcf} --psam ${sample_information_file} --output-chr chrM --make-bed --split-par ${plink_par} --out ${genotype}
    """

    stub:
    """
    touch ${genotype}.log
    touch ${genotype}.bed
    touch ${genotype}.bim
    touch ${genotype}.fam
    """
}