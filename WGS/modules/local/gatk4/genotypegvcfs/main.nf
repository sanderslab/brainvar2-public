process GATK4_GENOTYPEGVCFS {
    tag "$meta.contig"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container 'quay.io/biocontainers/gatk4:4.2.6.1--hdfd78af_0@sha256:2036bc919ce30eda27195ecd2639d6d556d4af03c04a66fcd5628211b2cd13db'

    input:
    tuple val(meta), path(genomicsdb)
    path fasta
    path fai
    path dict
    path dbsnp
    path dbsnp_index

    output:
    tuple val(meta), path('*.raw.vcf.gz'), path('*.raw.vcf.gz.tbi'), emit: vcf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.contig
    """
    gatk --java-options "-XX:MaxRAMPercentage=90" GenotypeGVCFs \\
        --variant gendb://${genomicsdb} \\
        --output ${prefix}.raw.vcf.gz \\
        --reference ${fasta} \\
        --dbsnp ${dbsnp} \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        --intervals ${meta.contig} \\
        --tmp-dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.contig
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' | gzip -c > ${prefix}.raw.vcf.gz
    printf 'stub index\n' > ${prefix}.raw.vcf.gz.tbi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: 4.2.6.1
    END_VERSIONS
    """
}
