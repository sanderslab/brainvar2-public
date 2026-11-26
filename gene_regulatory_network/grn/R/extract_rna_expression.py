#!/usr/bin/env python3

import sys
import pandas as pd
import anndata as ad
from scipy import sparse


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: python extract_rna_expression.py <input_h5ad> <gene_list.txt> <output.csv>",
            file=sys.stderr,
        )
        sys.exit(1)

    h5ad_file = sys.argv[1]
    gene_file = sys.argv[2]
    output_file = sys.argv[3]

    with open(gene_file, "r", encoding="utf-8") as f:
        genes = [line.strip() for line in f if line.strip()]

    if not genes:
        print("ERROR: No genes provided")
        sys.exit(1)

    adata = ad.read_h5ad(h5ad_file)

    # Case-insensitive gene symbol matching while preserving requested TF names.
    # This keeps motif TF labels (e.g. "Sox5") aligned with X matrix columns even
    # when h5ad stores symbols in a different case (e.g. "SOX5").
    var_names = pd.Index(adata.var_names.astype(str))
    var_upper = var_names.str.upper()

    upper_to_var = {}
    for i, g_up in enumerate(var_upper):
        if g_up not in upper_to_var:
            upper_to_var[g_up] = var_names[i]

    matched_pairs = []
    for req in genes:
        req_up = str(req).upper()
        if req_up in upper_to_var:
            matched_pairs.append((req, upper_to_var[req_up]))

    if len(matched_pairs) == 0:
        print(f"ERROR: No genes matched requested list ({len(genes)} genes)")
        sys.exit(1)

    # Use unique actual genes for efficient slicing, then remap to requested labels.
    unique_actual = []
    seen = set()
    for _, actual in matched_pairs:
        if actual not in seen:
            seen.add(actual)
            unique_actual.append(actual)

    adata_sub = adata[:, unique_actual]
    x = adata_sub.X
    if sparse.issparse(x):
        x = x.toarray()

    df_actual = pd.DataFrame(x, index=adata_sub.obs_names, columns=adata_sub.var_names)
    out_data = {}
    for req, actual in matched_pairs:
        out_data[req] = df_actual[actual].to_numpy()
    df = pd.DataFrame(out_data, index=adata_sub.obs_names)
    df.to_csv(output_file)

    print(f"Loaded h5ad: {h5ad_file}")
    print(f"Requested genes: {len(genes)}")
    print(f"Matched genes: {len(matched_pairs)}")
    print(f"Output shape: {df.shape[0]} rows x {df.shape[1]} cols")
    print(f"Saved: {output_file}")


if __name__ == "__main__":
    main()
