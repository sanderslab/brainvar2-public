#!/bin/bash
snakemake --use-conda --cores 1 --directory /workflow .snakemake/storage/s3/brainvar2/bulk_RNA_counts/htseq_count/$1/$2/htseq_counts.csv
