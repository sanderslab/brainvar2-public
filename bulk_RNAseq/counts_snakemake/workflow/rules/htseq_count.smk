import pandas as pd
paths = pd.DataFrame.from_records(
    (
        {"experiment": "BrainVar1", "path": "BrainVar2/Bulk_RNA/BAMS/BrainVar1"},
        {"experiment": "BrainVar2", "path": "BrainVar2/Bulk_RNA/BAMS"},
    ))

rule htseq_count:
    input:
        Aligned_sortedByCoord = collect("/data/sestan-lab-share/{item.path}/{{sample}}.Aligned.sortedByCoord.out.bam", item=lookup(query="experiment == '{experiment}'", within=paths))
    output:
        storage.s3(f"s3://{config['output_bucket']}/htseq_count/{{experiment}}/{{sample}}/htseq_counts.csv")
    log:
        storage.s3(f"s3://{config['output_bucket']}/htseq_count/{{experiment}}/{{sample}}/htseq_count.log")
    conda:
        "htseq"
    params:
        annotation_gtf=f"{config["annotation_gtf"]}",
        reverse_flag=branch(evaluate("{experiment} == 'BrainVar1'"), then="reverse", otherwise="yes")
    shell:
        "htseq-count -f bam -s {params.reverse_flag} -r pos -m union {input.Aligned_sortedByCoord} {params.annotation_gtf} -c {output} > {log} 2>&1"
