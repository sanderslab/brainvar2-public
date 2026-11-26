#!/usr/bin/env python3
import sys
import pandas as pd
df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
df = df.T.fillna(df.median(axis=1)).T
sorted=df.sort_index(axis=1)
sorted.to_csv(sys.argv[2], sep='\t', header=True, index=True)
