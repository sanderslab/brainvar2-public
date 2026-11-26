#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BRAINVAR2_WGS } from './workflows/brainvar2_wgs'

workflow {
    BRAINVAR2_WGS()
}
