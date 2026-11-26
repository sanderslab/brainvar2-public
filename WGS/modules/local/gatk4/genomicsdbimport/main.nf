process GATK4_GENOMICSDBIMPORT {
    tag "$meta.contig"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db'

    input:
    tuple val(meta), val(samples), path(gvcfs), path(gvcf_indexes)

    output:
    tuple val(meta), path('*.genomicsdb'), emit: genomicsdb
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.contig
    def args = task.ext.args ?: ''
    def sampleMap = samples.withIndex().collect { sample, index ->
        "${sample}\t${gvcfs[index]}"
    }.join('\n')
    """
cat > sample_map.tsv <<'END_SAMPLE_MAP'
${sampleMap}
END_SAMPLE_MAP

    gatk --java-options "-XX:MaxRAMPercentage=90" GenomicsDBImport \\
        --sample-name-map sample_map.tsv \\
        --genomicsdb-workspace-path ${prefix}.genomicsdb \\
        --intervals ${meta.contig} \\
        --tmp-dir . \\
        ${args}

cat > versions.yml <<'END_VERSIONS'
"${task.process}":
    gatk4: 4.2.6.1
END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.contig
    """
    mkdir -p ${prefix}.genomicsdb
    touch ${prefix}.genomicsdb/callset.json
cat > versions.yml <<'END_VERSIONS'
"${task.process}":
    gatk4: 4.2.6.1
END_VERSIONS
    """
}
