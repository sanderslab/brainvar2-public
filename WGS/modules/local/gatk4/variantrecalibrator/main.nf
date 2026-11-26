process GATK4_VARIANTRECALIBRATOR {
    tag "$meta.id:$mode"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db'

    input:
    tuple val(meta), path(vcf), path(vcf_index)
    path resources
    path resource_indexes
    val resource_specs
    path fasta
    path fai
    path dict
    val mode
    val annotations
    val tranches
    val max_gaussians

    output:
    tuple val(meta), path('*.recal'), path('*.recal.idx'), path('*.tranches'), emit: model
    path '*plots.R', optional: true, emit: plots
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${mode.toLowerCase()}"
    def resourceArgs = resources.withIndex().collect { resource, index ->
        def spec = resource_specs[index]
        "--resource:${spec.label},known=${spec.known},training=${spec.training},truth=${spec.truth},prior=${spec.prior} ${resource}"
    }.join(' ')
    def annotationArgs = annotations.collect { "--an ${it}" }.join(' ')
    def trancheArgs = tranches.collect { "--tranche ${it}" }.join(' ')
    def gaussianArgs = max_gaussians ? "--max-gaussians ${max_gaussians}" : ''
    """
    gatk --java-options "-XX:MaxRAMPercentage=90" VariantRecalibrator \\
        --variant ${vcf} \\
        --output ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        --reference ${fasta} \\
        --mode ${mode} \\
        ${resourceArgs} \\
        ${annotationArgs} \\
        ${trancheArgs} \\
        ${gaussianArgs} \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${mode.toLowerCase()}"
    """
    touch ${prefix}.recal
    touch ${prefix}.recal.idx
    touch ${prefix}.tranches
    touch ${prefix}_plots.R
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}
