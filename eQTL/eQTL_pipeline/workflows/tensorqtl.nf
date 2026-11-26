include {
   Transpose
   RunTensorqtl
   SortBedFile
   FillNaNWithRowMedian
   RunTensorqtlNominal
} from '../modules/tensorqtl'

workflow TEST {
   Channel
      .fromPath( params.sample_file )
      .splitCsv( header: true, sep: ',' )
      .map { row -> tuple( row.sample_group, row.pc_name, row.cell_group ) }
      .set {AllRowsChannel}

   AllRowsChannel.view()

   AllRowsChannel
      .map { sample_group, pc_name, cell_group ->
         bed_file = file("${params.work_path}/${sample_group}/${cell_group}.bed")
         sorted_bed_file = "${cell_group}_sorted.bed"
         target_folder = "${sample_group}"
         tuple(sample_group, pc_name, cell_group, bed_file, sorted_bed_file, target_folder)
      }
      .set {SortBedFileChannel}

   SortBedFileChannel.view(item -> "SortBedFileChannel: $item")

   SortBedFileChannel
      .map { sample_group, pc_name, cell_group, bed_file, sorted_bed_file, target_folder ->
         covariate_file = file("${params.work_path}/${sample_group}/${pc_name}/${cell_group}.cov_v2.txt")
         target_folder = "${sample_group}/${pc_name}"
         target_file = "${cell_group}.cov_v2_fill.txt"
         tuple(covariate_file, target_folder, target_file)
      }
      .set {FillNaNWithRowMedianChannel}

   FillNaNWithRowMedianChannel.view(item -> "FillNaNWithRowMedianChannel: $item")
}

workflow Prepare {

   Channel
      .fromPath( params.sample_file )
      .splitCsv( header: true, sep: ',' )
      .map { row -> tuple( row.sample_group, row.pc_name, row.cell_group, row.genotype ) }
      .set {AllRowsChannel}

   AllRowsChannel.view()

   AllRowsChannel
      .map { sample_group, pc_name, cell_group, genotype ->
         cov_file = file("${params.work_path}/01-raw/${sample_group}/${pc_name}/${cell_group}.cov.txt")
         transposed_file = "${cell_group}.cov.txt"
         target_folder = "inputs/${sample_group}/${pc_name}"
         tuple(sample_group, pc_name, cell_group, cov_file, transposed_file, target_folder)
      }
      .set {TransposeChannel}

   TransposeChannel.view(item -> "TransposeChannel: $item")
   TransposeOutput = Transpose(TransposeChannel)
}

workflow QTL {

   Channel
      .fromPath( params.sample_file )
      .splitCsv( header: true, sep: ',' )
      .map { row -> tuple( row.sample_group, row.pc_name, row.cell_group, row.genotype ) }
      .set {AllRowsChannel}

   AllRowsChannel.view()

   AllRowsChannel
      .map { sample_group, pc_name, cell_group, genotype ->
         bed_file = file("${params.work_path}/inputs/${sample_group}/${cell_group}.bed")
         sorted_bed_file = "${cell_group}_sorted.bed"
         target_folder = "inputs/${sample_group}"
         tuple(sample_group, pc_name, cell_group, genotype, bed_file, sorted_bed_file, target_folder)
      }
      .set {SortBedFileChannel}

   SortBedFileChannel.view(item -> "SortBedFileChannel: $item")
   SortBedFileOutput = SortBedFile(SortBedFileChannel)
   
   SortBedFileOutput.all_of_it
      .map { sample_group, pc_name, cell_group, genotype ->
         covariate_file = file("${params.work_path}/inputs/${sample_group}/${pc_name}/${cell_group}.cov.txt")
         target_folder = "inputs/${sample_group}/${pc_name}"
         target_file = "${cell_group}.cov_fill.txt"
         tuple(sample_group, pc_name, cell_group, genotype, covariate_file, target_folder, target_file)
      }
      .set {FillNaNWithRowMedianChannel}

   FillNaNWithRowMedianChannel.view(item -> "FillNaNWithRowMedianChannel: $item")
   FillNaNWithRowMedianOutput = FillNaNWithRowMedian (FillNaNWithRowMedianChannel)

   FillNaNWithRowMedianOutput.all_of_it
      .map { sample_group, pc_name, cell_group, genotype ->
         plink_bed = "${params.work_path}/plink_vcf/${genotype}.bed"
         plink_bim = "${params.work_path}/plink_vcf/${genotype}.bim"
         plink_fam = "${params.work_path}/plink_vcf/${genotype}.fam"
         covariates = file("${params.work_path}/inputs/${sample_group}/${pc_name}/${cell_group}.cov_fill.txt")
         phenotypes = file("${params.work_path}/inputs/${sample_group}/${cell_group}_sorted.bed")
         target_folder = "tensorqtl/${sample_group}/${pc_name}/${cell_group}"
         file_name = cell_group
         tuple(target_folder, file_name, genotype, plink_bed, plink_bim, plink_fam, covariates, phenotypes)
      }
      .set {TensorqtlChannel}

   TensorqtlOutput = RunTensorqtl (TensorqtlChannel)
   TensorqtlNominalOutput = RunTensorqtlNominal (TensorqtlChannel)
}

