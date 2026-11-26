from os import listdir
from os.path import isfile, join
def fastq_files(wildcards):
    mypath=f"/data/sestan-lab-share/BrainVar2/Bulk_RNA/FASTQ/Sample_{wildcards.sample}"
    onlyfiles = [join(mypath,f) for f in listdir(mypath) if isfile(join(mypath, f))]
    return onlyfiles

        # storage.s3("s3://sestan-lab-share/BrainVar2/Bulk_RNA/FASTQ/{sample}/{sample}_S??_L004_R1_001.fastq.gz"),
        # storage.s3("s3://sestan-lab-share/BrainVar2/Bulk_RNA/FASTQ/{sample}/{sample}_S??_L004_R2_001.fastq.gz")

rule salmon_counts:
    input:
        fastq_files
    output:
        quant_sf = storage.s3("s3://%s/salmon/BrainVar2/{sample}/quant.sf" % config["output_bucket"]),
        cmd_info = storage.s3("s3://%s/salmon/BrainVar2/{sample}/cmd_info.json" % config["output_bucket"]),
        lib_format_counts = storage.s3("s3://%s/salmon/BrainVar2/{sample}/lib_format_counts.json" % config["output_bucket"]),
        ambig_info = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/ambig_info.tsv" % config["output_bucket"]),
        expected_bias = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/expected_bias.gz" % config["output_bucket"]),
        fld = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/fld.gz" % config["output_bucket"]),
        meta_info = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/meta_info.json" % config["output_bucket"]),
        observed_bias = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/observed_bias.gz" % config["output_bucket"]),
        observed_bias_3p = storage.s3("s3://%s/salmon/BrainVar2/{sample}/aux_info/observed_bias_3p.gz" % config["output_bucket"]),
        flenDist = storage.s3("s3://%s/salmon/BrainVar2/{sample}/libParams/flenDist.txt" % config["output_bucket"]),
        log = storage.s3("s3://%s/salmon/BrainVar2/{sample}/logs/salmon_quant.log" % config["output_bucket"])
    log:
        storage.s3(f"s3://{config['output_bucket']}/salmon/BrainVar2/{{sample}}/salmon.log")
    params:
        output_folder=f".snakemake/storage/s3/{config["output_bucket"]}/salmon/BrainVar2/{{sample}}",
    conda:
        "salmon"
    shell:
        'salmon quant -i /data/salmon_index -l A -1 {input[0]} -2 {input[1]} --validateMappings -o {params.output_folder} > {log} 2>&1'
