
#!/bin/bash

# Automated GWAS-eQTL subsetting script - GENE-BASED FILTERING
# Filters eQTL files to genes identified in GWAS studies

# =============================================================================
# CONFIGURATION - MODIFY THESE PATHS
# =============================================================================

# Input directories
GWAS_DIR="/mnt/s3/eqtl/eQTL/GWAS_subset"
EQTL_DIR="/mnt/s3/eqtl/eQTL/munged"

# Output base directory
OUTPUT_BASE="/mnt/s3/eqtl/eQTL/subsetted_by_gene"

# eQTL file parameters
EQTL_DELIMITER=$'\t'
EQTL_EXT="tsv"

# Column names in eQTL files to match against (will try both)
GENE_SYMBOL_COL="GENESYMBOL"
ENSEMBL_COL="ENSEMBLID"

# =============================================================================
# GWAS NAME MAPPING
# =============================================================================

# Function to convert long GWAS filename to short abbreviation
get_gwas_abbreviation() {
    local filename=$1
    
    # Define mapping based on your R lists
    case "$filename" in
        *"WALK_metaGWAS_09032023"*) echo "AOW" ;;
        *"iPSYCH-PGC_ASD_Nov2017"*) echo "ASD" ;;
        *"PGC3_SCZ_wave3"*) echo "SCZ" ;;
        *"GCST90027158"*) echo "ALZ" ;;
        *"pgcAN2.2019-07"*) echo "AN" ;;
        *"GIANT_HEIGHT_YENGO_2022"*) echo "HEIGHT" ;;
        *"Hatoum2023AddictionEuropean"*) echo "SUD" ;;
        *"trans_ptsd_pcs_v4_aug3_2021"*) echo "PTSD" ;;
        *"daner_OCD_full_wo23andMe_190522"*) echo "OCD2025" ;;
        *"GCST90558311"*) echo "FTLD_TDP" ;;
        *"GCST90027163"*) echo "ALS" ;;
        *"TS_Oct2018"*) echo "TS" ;;
        *"PAU_trans_ancestry_Aug2023"*) echo "PAU" ;;
        *"bip2024_multianc_no23andMe"*) echo "BIP" ;;
        *"pgc-panic2019"*) echo "PANIC" ;;
        *"ADHD2022_iPSYCH_deCODE_PGC"*) echo "ADHD" ;;
        *"supplementary_data3_GWAS_meta"*) echo "DD" ;;
        *"GCST90271609"*) echo "GEPI" ;;
        *"GCST90271616"*) echo "CEPI" ;;
        *"EA4_additive_excl_23andMe"*) echo "EA" ;;
        *"ILAE_All_Epi_11.8.14"*) echo "EPI" ;;
        *"SavageJansen_2018_intelligence"*) echo "IQ" ;;
        *"pgc-mdd2025_no23andMe"*) echo "MDD" ;;
        *"discovery_metav3.0.meta"*) echo "MS" ;;
        *"ocd_aug2017"*) echo "OCD" ;;
        *"GCST90275127"*) echo "PD" ;;
        *) 
            # Fallback: extract a reasonable name from the filename
            local base=$(basename "$filename")
            base="${base#all_munged_}"
            base="${base#munged_}"
            base="${base%%INFOFilter*}"
            base="${base#70}"
            base="${base#no}"
            base="${base#_}"
            echo "$base"
            ;;
    esac
}

# =============================================================================
# FUNCTIONS
# =============================================================================

subset_eqtl_by_gene() {
    local gwas_file=$1
    local gwas_name=$2
    local eqtl_name=$3
    local eqtl_file=$4
    local output_dir=$5
    
    local output_file="${output_dir}/${eqtl_name}.${EQTL_EXT}"
    
    echo "  Processing: $eqtl_name with $gwas_name"
    
    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"
    
    # Use zcat or gzcat depending on system
    if command -v zcat &> /dev/null; then
        CAT_CMD="zcat"
    else
        CAT_CMD="gzcat"
    fi
    
    # Create a temporary file with the filtered results
    local temp_file="${output_dir}/temp_$$_${eqtl_name}.tsv"
    
    # Extract header and find column indices for gene identifiers
    local header=$($CAT_CMD "$eqtl_file" | head -n 1)
    
    # Find column numbers (awk uses 1-based indexing)
    local genesymbol_col=$(echo "$header" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="GENESYMBOL")print i}}')
    local ensembl_col=$(echo "$header" | awk -F'\t' '{for(i=1;i<=NF;i++){if($i=="ENSEMBLID")print i}}')
    
    if [ -z "$genesymbol_col" ] && [ -z "$ensembl_col" ]; then
        echo "    ERROR: Could not find GENESYMBOL or ENSEMBLID columns in eQTL file"
        return 1
    fi
    
    echo "    Using columns - GENESYMBOL: $genesymbol_col, ENSEMBLID: $ensembl_col"
    
    # Filter eQTL data by genes
    # This awk script:
    # 1. Loads all Ensembl IDs from GWAS file into memory
    # 2. Keeps the header line
    # 3. For each eQTL row, checks if GENESYMBOL or ENSEMBLID matches a GWAS gene
    $CAT_CMD "$eqtl_file" | awk -F'\t' -v gs_col="$genesymbol_col" -v ens_col="$ensembl_col" '
    BEGIN {OFS="\t"}
    NR==FNR {
        # Load Ensembl IDs from GWAS file (remove any whitespace and version numbers)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $1)
        # Remove version number from Ensembl ID (e.g., ENSG00000235098.5 -> ENSG00000235098)
        gene_id = $1
        sub(/\.[0-9]+$/, "", gene_id)
        genes[gene_id]=1
        genes[$1]=1  # Also keep original in case no version
        next
    }
    FNR==1 {
        # Print header
        print
        next
    }
    {
        # Remove version from eQTL Ensembl ID for matching
        if (gs_col != "") {
            gs_clean = $gs_col
            sub(/\.[0-9]+$/, "", gs_clean)
            if (gs_clean in genes || $gs_col in genes) {
                print
                next
            }
        }
        if (ens_col != "") {
            ens_clean = $ens_col
            sub(/\.[0-9]+$/, "", ens_clean)
            if (ens_clean in genes || $ens_col in genes) {
                print
                next
            }
        }
    }
    ' "$gwas_file" - > "$temp_file"
    
    # Check if output has more than just header
    local line_count=$(wc -l < "$temp_file")
    if [ "$line_count" -le 1 ]; then
        echo "    WARNING: No gene matches found for $eqtl_name with $gwas_name"
        rm "$temp_file"
        return
    fi
    
    # Move temp file to final location (header is already included)
    mv "$temp_file" "$output_file"
    
    # Compress output
    echo "    Compressing..."
    gzip -f "$output_file"
    
    local matches=$((line_count - 1))
    echo "    ✓ Created: ${output_file}.gz ($matches eQTL entries for matched genes)"
}

extract_eqtl_name() {
    local filename=$1
    local base=$(basename "$filename")
    base="${base%.tsv.gz}"
    
    # Remove common prefixes and suffixes
    base="${base#munged_eqtl_}"
    base="${base#munged_with_genes_eqtl_}"
    base="${base#gene_info_eqtl_}"
    base="${base%_pval_nominal_1e-5_chr}"
    
    echo "$base"
}

# =============================================================================
# MAIN PROCESSING
# =============================================================================

echo "=========================================="
echo "GWAS-eQTL Gene-Based Subsetting Script"
echo "=========================================="
echo ""

# Find all GWAS gene list files
echo "Searching for GWAS gene list files in: $GWAS_DIR"
GWAS_FILES=($(find "$GWAS_DIR" -type f -name "all_munged*_genes_in_2Mb_intervals.txt" | sort))

if [ ${#GWAS_FILES[@]} -eq 0 ]; then
    echo "ERROR: No GWAS gene list files found in $GWAS_DIR"
    echo "Looking for files matching pattern: all_munged*_genes_in_2Mb_intervals.txt"
    exit 1
fi

echo "Found ${#GWAS_FILES[@]} GWAS file(s)"

# Show first few genes from first GWAS file as verification
echo ""
echo "Sample genes from first GWAS file:"
head -5 "${GWAS_FILES[0]}" | sed 's/^/  /'
echo ""

# Find all eQTL files
echo "Searching for eQTL files in: $EQTL_DIR"
EQTL_FILES=($(find "$EQTL_DIR" -type f \( -name "munged_eqtl_*.tsv.gz" -o -name "munged_with_genes_eqtl_*.tsv.gz" \) | sort))

if [ ${#EQTL_FILES[@]} -eq 0 ]; then
    echo "ERROR: No eQTL files found in $EQTL_DIR"
    exit 1
fi

echo "Found ${#EQTL_FILES[@]} eQTL file(s):"
for eqtl_file in "${EQTL_FILES[@]}"; do
    eqtl_name=$(extract_eqtl_name "$eqtl_file")
    echo "  - $eqtl_name"
done
echo ""

# Calculate total combinations
total_combinations=$((${#GWAS_FILES[@]} * ${#EQTL_FILES[@]}))
current=0

echo "Will process $total_combinations combinations (${#GWAS_FILES[@]} GWAS × ${#EQTL_FILES[@]} eQTL)"
echo "Matching strategy: Filter eQTLs by GENESYMBOL or ENSEMBLID from GWAS gene lists"
echo "=========================================="
echo ""

# Process each GWAS file
for gwas_file in "${GWAS_FILES[@]}"; do
    # Extract GWAS name from filename using abbreviation mapping
    gwas_basename=$(basename "$gwas_file")
    gwas_name=$(get_gwas_abbreviation "$gwas_basename")
    gwas_full_name="${gwas_basename#all_munged_}"
    gwas_full_name="${gwas_full_name%.txt_genes_in_2Mb_intervals.txt}"
    
    echo "=========================================="
    echo "Processing GWAS: $gwas_name ($gwas_full_name)"
    echo "=========================================="
    
    # Count genes in this GWAS file
    gene_count=$(wc -l < "$gwas_file")
    echo "Gene count: $gene_count"
    echo ""
    
    # Create output directory for this GWAS
    output_dir="${OUTPUT_BASE}/${gwas_name}"
    
    # Process each eQTL file
    for eqtl_file in "${EQTL_FILES[@]}"; do
        current=$((current + 1))
        eqtl_name=$(extract_eqtl_name "$eqtl_file")
        
        echo "[$current/$total_combinations]"
        subset_eqtl_by_gene "$gwas_file" "$gwas_name" "$eqtl_name" "$eqtl_file" "$output_dir"
    done
    
    echo ""
done

echo "=========================================="
echo "All processing complete!"
echo "=========================================="
echo ""
echo "Output files are organized in: $OUTPUT_BASE"
echo "Structure: $OUTPUT_BASE/[GWAS_NAME]/[EQTL_NAME].tsv.gz"
echo ""
echo "Summary:"
echo "  - Processed ${#GWAS_FILES[@]} GWAS datasets"
echo "  - Processed ${#EQTL_FILES[@]} eQTL datasets"
echo "  - Filtered eQTLs to keep only genes identified in GWAS"
echo "  - Created up to $total_combinations output files"
