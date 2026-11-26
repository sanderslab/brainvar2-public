include { PARABRICKS_FQ2BAM } from '../modules/local/parabricks/fq2bam/main'
include { PARABRICKS_HAPLOTYPECALLER } from '../modules/local/parabricks/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT } from '../modules/local/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS } from '../modules/local/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_RAW } from '../modules/local/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as GATK4_MERGEVCFS_FINAL } from '../modules/local/gatk4/mergevcfs/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP } from '../modules/local/gatk4/variantrecalibrator/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL } from '../modules/local/gatk4/variantrecalibrator/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL } from '../modules/local/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP } from '../modules/local/gatk4/applyvqsr/main'

def validateKnownParams(parameters) {
    def known = [
        'stage', 'input', 'outdir', 'cohort', 'contigs',
        'parabricks_refdir', 'fasta_basename', 'known_sites',
        'known_sites_indexes', 'fasta', 'fai', 'dict', 'dbsnp',
        'dbsnp_index', 'hapmap', 'hapmap_index', 'omni', 'omni_index',
        'one_thousand_g', 'one_thousand_g_index', 'axiom', 'axiom_index',
        'mills', 'mills_index'
    ] as Set
    def unknown = parameters.keySet().findAll { !known.contains(it) }.sort()
    if (unknown) {
        error "Unknown workflow parameter(s): ${unknown.collect { '--' + it }.join(', ')}"
    }
}

def requireParam(value, name) {
    if (value == null || value.toString().trim() == '') {
        error "Missing required parameter: --${name}"
    }
    return value
}

def checkedFile(value, name) {
    requireParam(value, name)
    return file(value, checkIfExists: true)
}

def parseContigs(value) {
    def contigs = requireParam(value, 'contigs').toString().split(',')
        .collect { it.trim() }
        .findAll { it }
    if (!contigs || contigs.size() != contigs.toSet().size()) {
        error '--contigs must contain a non-empty, unique comma-separated list'
    }
    return contigs
}

def validateParabricksRow(row) {
    ['sample', 'fastq_1', 'fastq_2'].each { column ->
        if (!row.containsKey(column) || !row[column]?.trim()) {
            error "Parabricks samplesheet requires a non-empty '${column}' column"
        }
    }
    if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample identifier '${row.sample}'"
    }
    return tuple(
        [id: row.sample],
        [file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true)]
    )
}

def validateParabricksRows(rows) {
    if (!rows) {
        error 'The Parabricks samplesheet contains no data rows'
    }

    def counts = rows.countBy { meta, reads ->
        "${meta.id}\t${reads[0]}\t${reads[1]}"
    }
    def duplicates = counts.findAll { key, count -> count != 1 }.keySet().sort()
    if (duplicates) {
        error "Parabricks samplesheet contains duplicate sample-FASTQ records: ${duplicates.take(10).join(',')}"
    }
    return rows
}

def validateJointRows(rows, expectedContigs) {
    if (!rows) {
        error 'The joint-calling samplesheet contains no data rows'
    }

    def samples = rows.collect { it[0].id }.toSet().sort()
    def observedContigs = rows.collect { it[0].contig }.toSet()
    def unexpected = observedContigs - expectedContigs.toSet()
    if (unexpected) {
        error "Joint samplesheet contains contigs not requested by --contigs: ${unexpected.sort().join(',')}"
    }

    def counts = rows.countBy { "${it[0].id}\t${it[0].contig}" }
    def duplicates = counts.findAll { key, count -> count != 1 }.keySet().sort()
    if (duplicates) {
        error "Joint samplesheet has duplicate sample-contig records: ${duplicates.take(10).join(',')}"
    }

    def missing = []
    samples.each { sample ->
        expectedContigs.each { contig ->
            if (!counts.containsKey("${sample}\t${contig}")) {
                missing << "${sample}:${contig}"
            }
        }
    }
    if (missing) {
        error "Joint samplesheet is incomplete; missing records include: ${missing.take(10).join(',')}"
    }
    return rows
}

def jointRow(row) {
    ['sample', 'contig', 'gvcf', 'gvcf_index'].each { column ->
        if (!row.containsKey(column) || !row[column]?.trim()) {
            error "Joint samplesheet requires a non-empty '${column}' column"
        }
    }
    if (!(row.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
        error "Invalid sample identifier '${row.sample}'"
    }
    return tuple(
        [id: row.sample, contig: row.contig],
        file(row.gvcf, checkIfExists: true),
        file(row.gvcf_index, checkIfExists: true)
    )
}

def orderedCalls(rows, contigs, cohort) {
    def rank = contigs.withIndex().collectEntries { contig, index -> [(contig): index] }
    def ordered = rows.sort { left, right -> rank[left[0].contig] <=> rank[right[0].contig] }
    return tuple([id: cohort], ordered.collect { it[1] }, ordered.collect { it[2] })
}

// Contract: params-defined stage inputs; emits final cohort VCF or per-sample gVCFs plus version records.
workflow BRAINVAR2_WGS {
    main:
    validateKnownParams(params)
    requireParam(params.stage, 'stage')
    requireParam(params.input, 'input')

    ch_versions = Channel.empty()
    ch_final = Channel.empty()
    ch_parabricks_gvcf = Channel.empty()

    if (params.stage == 'parabricks') {
        def referenceDirectory = checkedFile(params.parabricks_refdir, 'parabricks_refdir')
        def knownSitePaths = requireParam(params.known_sites, 'known_sites').toString().split(',')
            .collect { checkedFile(it.trim(), 'known_sites') }
        def knownSiteIndexPaths = requireParam(params.known_sites_indexes, 'known_sites_indexes').toString().split(',')
            .collect { checkedFile(it.trim(), 'known_sites_indexes') }
        if (knownSitePaths.size() != knownSiteIndexPaths.size()) {
            error '--known_sites and --known_sites_indexes must contain the same number of files'
        }


        ch_reads = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row -> validateParabricksRow(row) }
            .collect(flat: false)
            .map { rows -> validateParabricksRows(rows) }
            .flatMap { rows -> rows }
            .map { meta, reads -> tuple(meta.id, reads) }
            .groupTuple()
            .map { sample, readPairs -> tuple([id: sample], readPairs.flatten()) }

        PARABRICKS_FQ2BAM(
            ch_reads,
            Channel.value(referenceDirectory),
            Channel.value(params.fasta_basename),
            Channel.value(knownSitePaths),
            Channel.value(knownSiteIndexPaths)
        )

        ch_haplotypecaller_input = PARABRICKS_FQ2BAM.out.cram
            .map { meta, cram, crai, recal, metrics -> tuple(meta, cram, crai, recal) }

        PARABRICKS_HAPLOTYPECALLER(
            ch_haplotypecaller_input,
            Channel.value(referenceDirectory),
            Channel.value(params.fasta_basename)
        )

        ch_parabricks_gvcf = PARABRICKS_HAPLOTYPECALLER.out.gvcf
        ch_versions = ch_versions
            .mix(PARABRICKS_FQ2BAM.out.versions.first())
            .mix(PARABRICKS_HAPLOTYPECALLER.out.versions.first())
    }
    else if (params.stage == 'joint') {
        def contigs = parseContigs(params.contigs)
        def cohort = requireParam(params.cohort, 'cohort')

        def fasta = checkedFile(params.fasta, 'fasta')
        def fai = checkedFile(params.fai, 'fai')
        def dict = checkedFile(params.dict, 'dict')
        def dbsnp = checkedFile(params.dbsnp, 'dbsnp')
        def dbsnpIndex = checkedFile(params.dbsnp_index, 'dbsnp_index')

        def snpResources = [
            checkedFile(params.hapmap, 'hapmap'),
            checkedFile(params.omni, 'omni'),
            checkedFile(params.one_thousand_g, 'one_thousand_g'),
            dbsnp
        ]
        def snpResourceIndexes = [
            checkedFile(params.hapmap_index, 'hapmap_index'),
            checkedFile(params.omni_index, 'omni_index'),
            checkedFile(params.one_thousand_g_index, 'one_thousand_g_index'),
            dbsnpIndex
        ]
        def snpSpecs = [
            [label: 'hapmap', known: false, training: true, truth: true, prior: 15.0],
            [label: 'omni', known: false, training: true, truth: true, prior: 12.0],
            [label: '1000G', known: false, training: true, truth: false, prior: 10.0],
            [label: 'dbsnp', known: true, training: false, truth: false, prior: 7.0]
        ]

        def indelResources = [
            checkedFile(params.axiom, 'axiom'),
            checkedFile(params.mills, 'mills'),
            dbsnp
        ]
        def indelResourceIndexes = [
            checkedFile(params.axiom_index, 'axiom_index'),
            checkedFile(params.mills_index, 'mills_index'),
            dbsnpIndex
        ]
        def indelSpecs = [
            [label: 'axiomPoly', known: false, training: true, truth: false, prior: 10.0],
            [label: 'mills', known: false, training: true, truth: true, prior: 12.0],
            [label: 'dbsnp', known: true, training: false, truth: false, prior: 2.0]
        ]

        ch_joint_rows = Channel
            .fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            // Cross-row validation requires the complete long-form manifest matrix.
            .map { row -> jointRow(row) }
            .collect(flat: false)
            .map { rows -> validateJointRows(rows, contigs) }
            .flatMap { rows -> rows }

        // Key: contig. Shape after grouping: [contig, samples[], gvcfs[], indexes[]].
        ch_by_contig = ch_joint_rows
            .map { meta, gvcf, index -> tuple(meta.contig, meta.id, gvcf, index) }
            .groupTuple(by: 0)
            .map { contig, samples, gvcfs, indexes ->
                tuple([id: contig, contig: contig], samples, gvcfs, indexes)
            }

        GATK4_GENOMICSDBIMPORT(ch_by_contig)

        GATK4_GENOTYPEGVCFS(
            GATK4_GENOMICSDBIMPORT.out.genomicsdb,
            Channel.value(fasta),
            Channel.value(fai),
            Channel.value(dict),
            Channel.value(dbsnp),
            Channel.value(dbsnpIndex)
        )

            // VQSR training requires the complete raw cohort callset across requested contigs.
        ch_raw_merged_input = GATK4_GENOTYPEGVCFS.out.vcf
            .collect(flat: false)
            .map { rows -> orderedCalls(rows, contigs, cohort) }

        GATK4_MERGEVCFS_RAW(ch_raw_merged_input, Channel.value(dict))

        GATK4_VARIANTRECALIBRATOR_SNP(
            GATK4_MERGEVCFS_RAW.out.vcf,
            Channel.value(snpResources),
            Channel.value(snpResourceIndexes),
            Channel.value(snpSpecs),
            Channel.value(fasta),
            Channel.value(fai),
            Channel.value(dict),
            Channel.value('SNP'),
            Channel.value(['QD', 'MQ', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR', 'DP']),
            Channel.value([100.0, 99.8, 99.6, 99.4, 99.2, 99.0, 95.0, 90.0]),
            Channel.value(0)
        )

        GATK4_VARIANTRECALIBRATOR_INDEL(
            GATK4_MERGEVCFS_RAW.out.vcf,
            Channel.value(indelResources),
            Channel.value(indelResourceIndexes),
            Channel.value(indelSpecs),
            Channel.value(fasta),
            Channel.value(fai),
            Channel.value(dict),
            Channel.value('INDEL'),
            Channel.value(['QD', 'MQRankSum', 'ReadPosRankSum', 'FS', 'SOR', 'DP']),
            Channel.value([100.0, 99.0, 95.0, 92.0, 90.0]),
            Channel.value(4)
        )

        ch_apply_indel = GATK4_GENOTYPEGVCFS.out.vcf
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.model)
            .map { callMeta, vcf, tbi, modelMeta, recal, recalIndex, tranches ->
                tuple(callMeta, vcf, tbi, recal, recalIndex, tranches)
            }

        GATK4_APPLYVQSR_INDEL(
            ch_apply_indel,
            Channel.value(fasta),
            Channel.value(fai),
            Channel.value(dict),
            Channel.value('INDEL'),
            Channel.value(99.0)
        )

        ch_apply_snp = GATK4_APPLYVQSR_INDEL.out.vcf
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.model)
            .map { callMeta, vcf, tbi, modelMeta, recal, recalIndex, tranches ->
                tuple(callMeta, vcf, tbi, recal, recalIndex, tranches)
            }

        GATK4_APPLYVQSR_SNP(
            ch_apply_snp,
            Channel.value(fasta),
            Channel.value(fai),
            Channel.value(dict),
            Channel.value('SNP'),
            Channel.value(99.8)
        )

        ch_final_merge_input = GATK4_APPLYVQSR_SNP.out.vcf
            // Final assembly requires all recalibrated contigs in explicit reference order.
            .collect(flat: false)
            .map { rows -> orderedCalls(rows, contigs, "${cohort}.jointcalled.recalibrated.GRCh38") }

        GATK4_MERGEVCFS_FINAL(ch_final_merge_input, Channel.value(dict))

        ch_final = GATK4_MERGEVCFS_FINAL.out.vcf
        ch_versions = ch_versions
            .mix(GATK4_GENOMICSDBIMPORT.out.versions.first())
            .mix(GATK4_GENOTYPEGVCFS.out.versions.first())
            .mix(GATK4_MERGEVCFS_RAW.out.versions)
            .mix(GATK4_VARIANTRECALIBRATOR_SNP.out.versions)
            .mix(GATK4_VARIANTRECALIBRATOR_INDEL.out.versions)
            .mix(GATK4_APPLYVQSR_INDEL.out.versions.first())
            .mix(GATK4_APPLYVQSR_SNP.out.versions.first())
            .mix(GATK4_MERGEVCFS_FINAL.out.versions)
    }
    else {
        error "Unsupported --stage '${params.stage}'. Choose 'parabricks' or 'joint'."
    }

    ch_versions
        .collectFile(
            name: 'software_versions.yml',
            newLine: true,
            storeDir: "${params.outdir}/pipeline_info"
        )

    emit:
    final_vcf = ch_final
    parabricks_gvcf = ch_parabricks_gvcf
    versions = ch_versions
}
