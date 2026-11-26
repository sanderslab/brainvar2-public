process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_gpu'

    container 'nvcr.io/nvidia/clara/clara-parabricks:4.0.1-1@sha256:1aed52b6730a213b8bf7a1aa6998668ad26f421fb49230c04e4340c6f68c4d9e'

    input:
    tuple val(meta), path(reads)
    path ref_dir
    val fasta_basename
    path known_sites
    path known_site_indexes

    output:
    tuple val(meta), path('*.cram'), path('*.cram.crai'), path('*.recal.table'), path('*.duplicate_metrics.txt'), emit: cram
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: meta.id
    def readPairs = reads.collate(2)
    if (!reads || reads.size() % 2 != 0) {
        error "PARABRICKS_FQ2BAM requires one or more complete FASTQ pairs"
    }
    def inFqArgs = readPairs.collect { pair ->
        "--in-fq ${pair[0]} ${pair[1]} " + '"$(read_group ' + pair[0] + ')"'
    }.join(' ')
    def knownSitesArgs = known_sites.collect { "--knownSites ${it}" }.join(' ')
    def numGpus = task.ext.num_gpus ?: 1
    // Preserve headroom for the Parabricks wrapper, compression, and BQSR
    // processes inside the task/container memory allocation.
    def memoryLimit = task.ext.parabricks_memory_gb ?: task.memory.giga / 2
    """
    read_group() {
        local read_one="\$1"
        local header flowcell lane
        set +o pipefail
        header=\$(gzip -cd "\$read_one" | head -n 1)
        set -o pipefail
        flowcell=\$(printf '%s\n' "\$header" | cut -d: -f3)
        lane=\$(printf '%s\n' "\$header" | cut -d: -f4)
        printf '%s' "@RG\\tID:\${flowcell}.${prefix}.\${lane}\\tPU:\${flowcell}.\${lane}\\tSM:${prefix}\\tLB:${prefix}\\tDS:${fasta_basename}\\tPL:illumina"
    }

    pbrun fq2bam \\
        --ref ${ref_dir}/${fasta_basename} \\
        ${inFqArgs} \\
        --out-duplicate-metrics ${prefix}.duplicate_metrics.txt \\
        ${knownSitesArgs} \\
        --out-recal-file ${prefix}.recal.table \\
        --bwa-options="-Y" \\
        --fix-mate \\
        --markdups-picard-version-2182 \\
        --optical-duplicate-pixel-distance 2500 \\
        --num-gpus ${numGpus} \\
        --memory-limit ${memoryLimit} \\
        --out-bam ${prefix}.deduped.sorted.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: 4.0.1-1
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch ${prefix}.deduped.sorted.cram
    touch ${prefix}.deduped.sorted.cram.crai
    touch ${prefix}.recal.table
    touch ${prefix}.duplicate_metrics.txt
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parabricks: 4.0.1-1
    END_VERSIONS
    """
}
