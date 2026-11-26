include { RunPlink } from '../modules/plink'

workflow PLINK {

   Channel
      .fromPath( params.genotypes_file )
      .splitCsv( header: true, sep: ',' )
      .map { row -> row.genotype }
      .set {AllRowsChannel}

   AllRowsChannel.view(item -> "AllRowsChannel: $item")

   AllRowsChannel
      .map { genotype ->
         vcf_file = "${params.genotype_location}/${genotype}.vcf.gz"
         sample_information_file = "${params.genotype_location}/${genotype}.samples.psam"
         plink_par = params.plink_par
         tuple(genotype, file(vcf_file), file(sample_information_file), plink_par)
      }
      .set {VcfFiles}

   VcfFiles.view(item -> "VcfFiles: $item")
   plink_output_ch = RunPlink (VcfFiles)

}
