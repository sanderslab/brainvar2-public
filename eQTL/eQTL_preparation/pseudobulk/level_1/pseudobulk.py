# Mount cell count matrices
goofys --region us-east-1 sestan-lab-share:BrainVar2/ data/

import anndata
import pandas as pd
import scanpy as sc
import decoupler as dc
import numpy as np

# Load the h5ad file
adata = anndata.read_h5ad('/home/ec2-user/data/Processed/level_1/BrainVar2_period_9_level_1.h5ad')

# Load the annotation/meta-data file into a DataFrame
annotation_file = "/home/ec2-user/data/Processed/level_1/BrainVar2_level1_cell_meta.csv"
cell_metadata = pd.read_csv(annotation_file)

# Merge the annotation DataFrame with the obs DataFrame of adata
merged_data = pd.merge(adata.obs, cell_metadata, left_on='Barcode', right_on='cell_id', how='left')

# Sort merged_data based on the index of adata.obs
merged_data = merged_data.set_index('Barcode')

# Update the 'cell_type' column in the obs DataFrame of adata with the merged cell type information
adata.obs['celltype'] = merged_data['celltype']
adata.obs['library'] = merged_data['library']

adata.X = np.round(adata.X).astype(int)
adata.layers['counts'] = adata.X

# Here, "mean" can be changed to "sum"
pdata = dc.get_pseudobulk(
    adata,
    sample_col='library',
    groups_col='celltype',
    layer='counts',
    mode='mean',
    min_cells=0,
    min_counts=0,
    min_smpls=0,
)

# Access the pseudobulk counts for each gene
pseudobulk_counts = pdata.X

# Convert the pseudobulk counts to a Pandas DataFrame
pseudobulk_counts_df = pd.DataFrame(data=pseudobulk_counts, columns=pdata.var_names, index=pdata.obs_names)

pseudobulk_counts_df.to_csv('/home/ec2-user/BrainVar2_period_9_level_1_pseudobulk_mean.txt', sep='\t')
