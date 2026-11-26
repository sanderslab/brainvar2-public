#!/bin/bash

# Path to the input file with sample information
INPUT_FILE="/home/ec2-user/file_lists/new_RNA_libraries.txt"

# Loop through the input file
while read -r R1 R2 SAMPLE_NAME; do
    OUTPUT_DIR=/home/ec2-user/STAR_alignment/"$SAMPLE_NAME"
    
    # Run STAR alignment
	STAR --genomeDir /home/ec2-user/STAR_genome_GTEx \
	     --runThreadN 14 \
	     --readFilesCommand zcat \
	     --readFilesIn $R1 $R2 \
	     --outFileNamePrefix $OUTPUT_DIR/
			 --quantMode GeneCounts \
			 --outSAMtype BAM SortedByCoordinate \

done < "$INPUT_FILE"
