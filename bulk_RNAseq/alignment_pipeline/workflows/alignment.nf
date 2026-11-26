include {
   RunStar
} from '../modules/star'

workflow Alignment {

   Channel
      .fromPath( params.sample_file )
      .splitCsv( header: true, sep: ',' )
      .map { row -> tuple( row.project, row.sample, row.file_path ) }
      .set {AllRowsChannel}

   AllRowsChannel.view()

   AllRowsChannel
      .map { project, sample, file_path ->
        def r1, r2
        if (project == "BrainVar2") {
            r1 = file("${params.s3_bucket}/${file_path}_R1_001.fastq.gz")
            r2 = file("${params.s3_bucket}/${file_path}_R2_001.fastq.gz")
        } 
        else {
            r1 = file("${params.s3_bucket}/${file_path}.R1.fastq.gz")
            r2 = file("${params.s3_bucket}/${file_path}.R2.fastq.gz")
        }
         star_index = file("${params.star_index}")
         gtf_file = file("${params.gtf_file}")
         target_folder = "${params.aligned_folder}/${project}/${sample}"
         tuple(sample, r1, r2, star_index, gtf_file, target_folder)
      }
      .set {AlignmentChannel}

   AlignmentChannel.view(item -> "AlignmentChannel: $item")
   AlignmentOutput = RunStar(AlignmentChannel)

}

