process PARABRICKS_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_gpu'

    container 'nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e'

    input:
    tuple val(meta), path(cram), path(crai), path(recal)
    path ref_dir
    val fasta_basename

    output:
    tuple val(meta), path('*.g.vcf.gz'), path('*.g.vcf.gz.tbi'), emit: gvcf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def numGpus = task.ext.num_gpus ?: 1
    """
    pbrun haplotypecaller \\
        --ref ${ref_dir}/${fasta_basename} \\
        --in-bam ${cram} \\
        --gvcf \\
        -G StandardAnnotation \\
        -G StandardHCAnnotation \\
        -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 \\
        -GQB 60 -GQB 70 -GQB 80 -GQB 90 \\
        --in-recal-file ${recal} \\
        --num-gpus ${numGpus} \\
        --out-variants ${prefix}.grch38.gatk.haplotypecaller.g.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: 4.0.1-1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > ${prefix}.grch38.gatk.haplotypecaller.g.vcf.gz
    printf 'stub index\n' > ${prefix}.grch38.gatk.haplotypecaller.g.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: 4.0.1-1
    END_VERSIONS
    """
}
