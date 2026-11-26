import pandas as pd
import gzip
from pyliftover import LiftOver
import re

# Initialize LiftOver converter
lo = LiftOver('hg38', 'hg19')

# File paths
files = [
    "/mnt/s3/eqtl/sQTL/munged/munged_sqtl_postnatal_1e-5.tsv.gz",
    "/mnt/s3/eqtl/sQTL/munged/munged_sqtl_prenatal_1e-5.tsv.gz"
]

def process_sqtl_liftover(infile, lo):
    """
    Process one sQTL file and perform liftover from hg38 to hg19/hg37
    """
    print(f"Processing: {infile}")
    
    # Read the gzipped TSV file
    dt = pd.read_csv(infile, sep='\t', compression='gzip')
    
    # Extract hg38 coordinates from PHENOTYPE_ID (format: chr:start:end:clu_id)
    coords = dt['PHENOTYPE_ID'].str.split(':', expand=True)
    
    # Add 'chr' prefix if not present and correct intervals
    dt['chr_hg38'] = coords[0].apply(lambda x: x if x.startswith('chr') else f'chr{x}')
    dt['start_hg38'] = pd.to_numeric(coords[1]).combine(pd.to_numeric(coords[2]), min)
    dt['end_hg38'] = pd.to_numeric(coords[1]).combine(pd.to_numeric(coords[2]), max)
    
    # Initialize hg37 columns
    dt['chr_hg37'] = None
    dt['start_hg37'] = None
    dt['end_hg37'] = None
    
    # Track unmapped coordinates and rows to append for multiple mappings
    unmapped_indices = []
    rows_to_append = []
    
    # Perform liftover for each coordinate
    for idx, row in dt.iterrows():
        chr_hg38 = row['chr_hg38']
        start_hg38 = int(row['start_hg38'])
        end_hg38 = int(row['end_hg38'])
        
        # LiftOver start and end positions (returns list of possible mappings)
        start_converted = lo.convert_coordinate(chr_hg38, start_hg38)
        end_converted = lo.convert_coordinate(chr_hg38, end_hg38)
        
        # Check if both positions mapped successfully
        if start_converted and end_converted:
            # Find all valid mappings where start and end map to same chromosome
            valid_mappings = []
            
            for start_result in start_converted:
                start_chr, start_pos = start_result[0], start_result[1]
                for end_result in end_converted:
                    end_chr, end_pos = end_result[0], end_result[1]
                    
                    # Only keep if both map to same chromosome
                    if start_chr == end_chr:
                        valid_mappings.append({
                            'chr': start_chr,
                            'start': min(start_pos, end_pos),
                            'end': max(start_pos, end_pos)
                        })
            
            if valid_mappings:
                # Add first mapping to original row
                dt.at[idx, 'chr_hg37'] = valid_mappings[0]['chr']
                dt.at[idx, 'start_hg37'] = valid_mappings[0]['start']
                dt.at[idx, 'end_hg37'] = valid_mappings[0]['end']
                
                # Create duplicate rows for additional mappings
                if len(valid_mappings) > 1:
                    for mapping in valid_mappings[1:]:
                        new_row = row.copy()
                        new_row['chr_hg37'] = mapping['chr']
                        new_row['start_hg37'] = mapping['start']
                        new_row['end_hg37'] = mapping['end']
                        rows_to_append.append(new_row)
            else:
                unmapped_indices.append(idx)
        else:
            unmapped_indices.append(idx)
    
    # Append all duplicate rows for multiple mappings
    if rows_to_append:
        dt = pd.concat([dt, pd.DataFrame(rows_to_append)], ignore_index=True)
        print(f"{len(rows_to_append)} additional rows created for multiple mappings")
    
    # Log and save unmapped loci
    if unmapped_indices:
        print(f"{len(unmapped_indices)} coordinates failed liftover")
        unmapped_file = infile.replace('.tsv.gz', '_unmapped.tsv')
        dt.iloc[unmapped_indices].to_csv(unmapped_file, sep='\t', index=False)
        print(f"Unmapped coordinates saved to: {unmapped_file}")
    
    # Save output
    outfile = infile.replace('.tsv.gz', '_hg37.tsv.gz')
    dt.to_csv(outfile, sep='\t', index=False, compression='gzip')
    print(f"Saved: {outfile}")
    
    return dt

# Process all files
results = []
for file in files:
    result = process_sqtl_liftover(file, lo)
    results.append(result)

print("All files processed successfully!")
