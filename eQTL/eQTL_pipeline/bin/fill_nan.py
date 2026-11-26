#!/usr/bin/env python3
import sys
import pandas as pd
import re

def is_text(value):
    value_str = str(value)
    pattern = re.compile(r'[a-zA-Z]')
    return bool(pattern.search(value_str))

def has_string_data(df, sample_size=10):
    string_columns = []
    for col in df.columns:
        sample = df[col].dropna().head(sample_size)
        if any(is_text(val) for val in sample):
            string_columns.append(col)
    return string_columns


df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)

dft = df.T
text_cols = has_string_data(dft)
cols = dft.columns.drop(text_cols)
dft[cols]=dft[cols].apply(pd.to_numeric, errors='coerce')
dft[cols]=dft[cols].fillna(dft[cols].median(axis=0))

df = dft[cols].T
sorted=df.sort_index(axis=1)
sorted.to_csv(sys.argv[2], sep='\t', header=True, index=True)
