process GATK4_MERGEVCFS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db'

    input:
    tuple val(meta), path(vcfs), path(vcf_indexes)
    path dict

    output:
    tuple val(meta), path('*.vcf.gz'), path('*.vcf.gz.tbi'), emit: vcf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def inputs = vcfs.collect { "--INPUT ${it}" }.join(' ')
    """
    gatk --java-options "-XX:MaxRAMPercentage=90" MergeVcfs \\
        ${inputs} \\
        --OUTPUT ${prefix}.vcf.gz \\
        --SEQUENCE_DICTIONARY ${dict} \\
        --CREATE_INDEX true \\
        --TMP_DIR .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > ${prefix}.vcf.gz
    printf 'stub index\n' > ${prefix}.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}
