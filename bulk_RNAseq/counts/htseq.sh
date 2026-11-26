
#!/bin/bash

# Iterate over every file ending in ".bam" in the current directory
for bam_file in *.bam; do
    # Extract the portion of the filename before "Aligned"
    output_file=$(basename "$bam_file" | sed 's/Aligned.*//' | sed 's/\./_/' | sed 's/_$//')"_counts.csv"

    # Run htseq-count command with the constructed output file name
    htseq-count -f bam -s yes -m union ~/bv2_bams/$bam_file ~/ref/gencode.v41.annotation.gtf -c ~/bv2_counts/$output_file
done

cd ~/bv1_bams/

# Iterate over every file ending in ".bam" in the current directory
for bam_file in *.bam; do
    # Extract the portion of the filename before "Aligned"
    output_file=$(basename "$bam_file" | sed 's/Aligned.*//' | sed 's/\./_/' | sed 's/_$//')"_counts.csv"

    # Run htseq-count command with the constructed output file name
    htseq-count -f bam -s yes -m union ~/bv1_bams/$bam_file ~/ref/gencode.v41.annotation.gtf -c ~/bv1_counts/$output_file
done
