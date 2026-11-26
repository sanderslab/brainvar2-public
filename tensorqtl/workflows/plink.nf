include { RUN_PLINK } from '../modules/plink'

workflow PLINK {

   Channel
      .value([params.genotype_vcf, params.sample_information_file, params.plink_par, params.genotype_plink_name])
      .set {vcf_channel}

   vcf_channel.view { tuple -> "Value: ${tuple[0]}, File: ${tuple[1]}" }
   plink_output_ch = RUN_PLINK (vcf_channel)

}
