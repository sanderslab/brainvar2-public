nextflow.enable.dsl=2

include { QTL; TEST; Prepare } from './workflows/tensorqtl'
include { PLINK } from './workflows/plink'

workflow {

    if (params.mode == "test") {
        TEST()
        } 
    else if (params.mode == "prepare") {
        Prepare()
        } 
    else if (params.mode == "qtl") {
        QTL()
        } 
    else if (params.mode == "plink") {
        PLINK()
        } 
    else {
        error "Invalid mode: ${params.mode}. Use 'qtl'"
    }

}
