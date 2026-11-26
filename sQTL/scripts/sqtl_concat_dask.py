#!/usr/bin/env python3
import pandas as pd
import glob
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import os
import numpy as np

def load_parquets(input_pattern):
    files = glob.glob(input_pattern)
    if not files:
        raise FileNotFoundError(f"No files match the pattern: {input_pattern}")
    print(f"Found {len(files)} files.")
    dfs = [pd.read_parquet(f) for f in files]
    return pd.concat(dfs, ignore_index=True)

def adjust_pvalues(df, pval_column='pval_nominal'):
    df[pval_column] = df[pval_column].fillna(1.0)
    df['pval_adj'] = multipletests(df[pval_column], method='fdr_bh')[1]
    return df

def create_qqplot(pvals, output_file):
    pvals = pvals.dropna()
    expected = -np.log10(np.linspace(1/len(pvals), 1, len(pvals)))
    observed = -np.log10(np.sort(pvals))
    plt.figure(figsize=(6,6))
    plt.scatter(expected, observed, s=10)
    plt.plot([0, max(expected)], [0, max(expected)], color='red', linestyle='--')
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    plt.title('QQ Plot')
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"QQ plot saved to {output_file}")

def create_pval_effect_plot(df, pval_column, effect_column, output_file):
    plt.figure(figsize=(8,6))
    plt.scatter(df[effect_column], -np.log10(df[pval_column].fillna(1.0)), s=5, alpha=0.5)
    plt.xlabel(effect_column)
    plt.ylabel(f"-log10({pval_column})")
    plt.title(f"{effect_column} vs {pval_column}")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()
    print(f"P-value vs effect size plot saved to {output_file}")

def process_dataset(dataset):
    print(f"Processing {dataset['name']} dataset")
    df = load_parquets(dataset['input_pattern'])

    # Filter invalid variant_ids
    df = df[~df['variant_id'].isin(['NA', '.'])]
    if 'ma_count' in df.columns:
        df = df[df['ma_count']>2]

    # Adjust p-values
    df = adjust_pvalues(df, pval_column=dataset['pval_column'])

    # Split phenotype_id
    split_cols = df['phenotype_id'].str.split(':', n=3, expand=True)
    split_cols.columns = ['chr', 'start', 'end', 'clu_info']
    df[['chr', 'start', 'end', 'clu_info']] = split_cols
    df['GeneSymbol'] = df['clu_info']
    df['GeneId'] = df['phenotype_id']

    # Create plots
    create_qqplot(df[dataset['pval_column']], dataset['qqplot_path'])
    create_pval_effect_plot(df, dataset['pval_column'], dataset['effect_column'], dataset['effect_plot_path'])

    # Filter by pval threshold
    if dataset['pval_threshold'] is not None:
        df = df[df[dataset['pval_column']] < dataset['pval_threshold']]
        print(f"{len(df)} rows pass {dataset['pval_column']} < {dataset['pval_threshold']}")

    # Save single parquet file
    os.makedirs(os.path.dirname(dataset['output_path']), exist_ok=True)
    df.to_parquet(dataset['output_path'], index=False)
    print(f"Saved filtered data to {dataset['output_path']}\n")

if __name__ == "__main__":
    datasets = [
        {
            'name': 'prenatal',
            'input_pattern': '/mnt/s3/eqtl/sQTL/tensorqtl/Common/Prenatal/PC8/Prenatal_250730/cis_nominal/*.parquet',
            'output_path': './results_sQTL/sQTL_prenatal_filtered.parquet',
            'pval_column': 'pval_nominal',
            'effect_column': 'slope',
            'pval_threshold': 0.05,
            'qqplot_path': './results_sQTL/qq_prenatal.png',
            'effect_plot_path': './results_sQTL/pval_effect_prenatal.png'
        },
        {
            'name': 'postnatal',
            'input_pattern': '/mnt/s3/eqtl/sQTL/tensorqtl/Common/Postnatal/PC8/Postnatal_250730/cis_nominal/*.parquet',
            'output_path': './results_sQTL/sQTL_postnatal_filtered.parquet',
            'pval_column': 'pval_nominal',
            'effect_column': 'slope',
            'pval_threshold': 0.05,
            'qqplot_path': './results_sQTL/qq_postnatal.png',
            'effect_plot_path': './results_sQTL/pval_effect_postnatal.png'
        }
    ]

    for dataset in datasets:
        process_dataset(dataset)
