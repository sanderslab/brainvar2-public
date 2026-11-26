process GATK4_APPLYVQSR {
    tag "$meta.contig:$mode"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db'

    input:
    tuple val(meta), path(vcf), path(vcf_index), path(recal), path(recal_index), path(tranches)
    path fasta
    path fai
    path dict
    val mode
    val truth_sensitivity

    output:
    tuple val(meta), path('*.vqsr.vcf.gz'), path('*.vqsr.vcf.gz.tbi'), emit: vcf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.contig}.${mode.toLowerCase()}"
    """
    gatk --java-options "-XX:MaxRAMPercentage=90" ApplyVQSR \\
        --variant ${vcf} \\
        --output ${prefix}.vqsr.vcf.gz \\
        --reference ${fasta} \\
        --tranches-file ${tranches} \\
        --recal-file ${recal} \\
        --mode ${mode} \\
        --truth-sensitivity-filter-level ${truth_sensitivity} \\
        --create-output-variant-index true \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.contig}.${mode.toLowerCase()}"
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > ${prefix}.vqsr.vcf.gz
    printf 'stub index\n' > ${prefix}.vqsr.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}
