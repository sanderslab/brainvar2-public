nextflow.enable.dsl=2

include { Alignment } from './workflows/alignment'

workflow {

    if (params.mode == "align") {
        Alignment()
        } 
    else {
        error "Invalid mode: ${params.mode}. Use 'align'"
    }

}
