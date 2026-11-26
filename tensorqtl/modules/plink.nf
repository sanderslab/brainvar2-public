process RUN_PLINK {
    label 'plink'
    containerOptions = '--env OUTPUT_DESTINATION=ox.paediatrics.sanders.batch'  
    publishDir "${params.work_path}/plink_vcf", mode: 'copy'

    input:
    tuple path(genotype_vcf), path(sample_information_file), val(plink_par), val(output_prefix)
   
    output:
    path "${output_prefix}.log", emit: plink_log
    path "${output_prefix}.bed", emit: genotype_bed
    path "${output_prefix}.bim", emit: variant_information
    path "${output_prefix}.fam", emit: sample_information

    script:
    """
    plink2 --vcf ${genotype_vcf} --psam ${sample_information_file} --output-chr chrM --make-bed --split-par ${plink_par} --out ${output_prefix}
    """

    stub:
    """
    touch ${output_prefix}.log
    touch ${output_prefix}.bed
    touch ${output_prefix}.bim
    touch ${output_prefix}.fam
    """
}