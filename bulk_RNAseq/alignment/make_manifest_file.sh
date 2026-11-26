find /home/ec2-user/fastqs -type f -name '*.fastq.gz' | sort > /home/ec2-user/file_lists/RNA_file_manifest.txt

#!/bin/bash

# Input file
input_file="/home/ec2-user/file_lists/RNA_file_manifest.txt"
# old_RNA_file_manifest.txt

# Output file
output_file="/home/ec2-user/file_lists/RNA_libraries.txt"
# old_RNA_libraries.txt

# Initialize the output file
echo -n > "$output_file"

# Read the input file line by line
while IFS= read -r r1_fastq; do
  # Read the next line (R2 fastq)
  read -r r2_fastq

  # Extract the library name from the R1 fastq file
  library_name=$(basename "$r1_fastq" | awk -F "Sample_" '{print $1}'| awk -F "_" '{print $1}')

  # Format the output line with tab-delimited fields
  output_line="$r1_fastq\t$r2_fastq\t$library_name"

  # Write the formatted entry to the output file
  echo -e "$output_line" >> "$output_file"
done < "$input_file"
