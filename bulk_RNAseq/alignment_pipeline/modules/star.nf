process RunStar {
    label 'star'
    publishDir "${params.work_path}", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2), path(star_index), path(gtf_file), val(target_folder)
   
    output:
    path "${target_folder}/*"

    script:
    def file_name_prefix = "${target_folder}/"
    """
    STAR --genomeDir ${star_index} \
        --runThreadN ${params.threads} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand gunzip -c \
        --outSAMattributes NH HI AS NM MD nM \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --sjdbGTFfile ${gtf_file} \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --sjdbScore 1 \
        --outFilterMismatchNmax 999 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --outFileNamePrefix ${file_name_prefix} \
        --limitSjdbInsertNsj 1000000

    """

    stub:
    """
    touch "${target_folder}_"
    """
}
