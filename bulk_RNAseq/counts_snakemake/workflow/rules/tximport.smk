rule tximport:
    input:
        rules.salmon_counts.output.quant_sf
    output:
        storage.s3(f"s3://{config['output_bucket']}/tximport/BrainVar2/{{sample}}/txi_counts.csv")
    log:
        storage.s3(f"s3://{config['output_bucket']}/tximport/BrainVar2/{{sample}}/tximport.log")
    script:
        "../scripts/tximport.R"