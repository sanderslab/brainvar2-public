#!/bin/bash

# April 2, 2025 by Narjes Rohani

# Directory containing the duplication rate text files calculated by 2_brainvar1_duplication_picard_calculation_batch

input_directory=".//bulkRNA_QC/code_picard/brainvar1/picard_dup_results_brainvar1"

# Output file to store the PERCENT_DUPLICATION values
output_file=".//bulkRNA_QC/code_picard/brainvar1/summary_percent_duplication_values.txt"

# Write the header to the output file
echo -e "Filename\tPERCENT_DUPLICATION" > "$output_file"

# Loop through each .txt file in the input directory
for file in "$input_directory"/*.txt; do
    # Check if the file contains the PERCENT_DUPLICATION keyword and extract the next line
    while read -r line; do
        # Look for the line with PERCENT_DUPLICATION
        if [[ "$line" =~ PERCENT_DUPLICATION ]]; then
            # Read the next line (where the duplication percentage is located)
            read -r next_line
            # Extract the PERCENT_DUPLICATION value from the next line
            percent_duplication=$(echo "$next_line" | awk '{print $9}')  # Assuming PERCENT_DUPLICATION is the 9th field
            
            # If a PERCENT_DUPLICATION value is found, write to the output file
            if [[ -n "$percent_duplication" ]]; then
                filename=$(basename "$file")
                echo -e "$filename\t$percent_duplication" >> "$output_file"
            fi
        fi
    done < "$file"
done

echo "Percent duplication values have been saved to $output_file"

