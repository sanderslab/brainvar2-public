import pandas as pd
import glob
import argparse
from statsmodels.stats.multitest import multipletests
import os

def load_and_concatenate_parquets(input_pattern):
    files = glob.glob(input_pattern)
    if not files:
        raise FileNotFoundError(f"No files match the pattern: {input_pattern}")
    print(f"Found {len(files)} files.")
    print(files)
    df_list = [pd.read_parquet(f) for f in files]
    return pd.concat(df_list, ignore_index=True)

def adjust_pvalues(df, pval_column='pval', method='fdr_bh'):
    if pval_column not in df.columns:
        raise KeyError(f"'{pval_column}' column not found in DataFrame.")
    adjusted_pvals = multipletests(df[pval_column], method=method)[1]
    df['pval_adj'] = adjusted_pvals
    return df

def main(input_pattern, output_path,pval_path, pval_column, fdr_threshold=None):
    # Load data
    df = load_and_concatenate_parquets(input_pattern)
    print('Concatenated data.')

    # Remove rows where variant_id is NA or "."

    df = df[~df['variant_id'].isin(['NA', '.'])]
    df = df[df['ma_samples']>2]

    print(f"{len(df)} rows after filtering invalid variant_ids.")

    # Fill missing p-values with 1
    df[pval_column] = df[pval_column].fillna(1.0)

    # Adjust p-values
    df = adjust_pvalues(df, pval_column=pval_column)
    pval_df = df[[pval_column,'pval_adj','slope']]
    print('Adjusted p-values calculated.')
    df = df[df[pval_column]<0.05]

    # Split phenotype_id into components
    # Format: chr:start:end:clu_ID_info
    df[['chr', 'start', 'end', 'clu_info']] = df['phenotype_id'].str.split(':', expand=True)
    df['GeneSymbol'] = df['clu_info']  # or extract differently if needed
    df['GeneId'] = df['phenotype_id']

    # Filter for significant results if threshold is set
    if fdr_threshold is not None:
        df = df[df['pval_adj'] < fdr_threshold]
        print(f"{len(df)} rows pass FDR < {fdr_threshold}")

    # Save result
    os.makedirs(os.path.dirname("./results_sQTL"), exist_ok=True)
    df.to_parquet(output_path, index=False)
    pval_df.to_csv(pval_path)
    print(f"Saved result to {output_path}")

if __name__ == "__main__":
    # Example: prenatal
    main(
        input_pattern='/mnt/s3/eqtl/sQTL/tensorqtl/Common/Prenatal/PC8/Prenatal_250730/cis_nominal/*.parquet',
        output_path='./results_sQTL/sQTL_prenatal_all_chr.parquet',
        pval_path = './results_sQTL/sQTL_prenatal_all_chr.csv',
        pval_column='pval_nominal',
        fdr_threshold=None
    )

    # Example: postnatal
    main(
        input_pattern='/mnt/s3/eqtl/sQTL/tensorqtl/Common/Postnatal/PC8/Postnatal_250730/cis_nominal/*.parquet',
        output_path='./results_sQTL/sQTL_postnatal_all_chr.parquet',
        pval_path = './results_sQTL/sQTL_posnatal_all_chr.csv',
        pval_column='pval_nominal',
        fdr_threshold=None
    )
