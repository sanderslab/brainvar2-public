#!/usr/bin/env python3
import sys
import pandas as pd

df = pd.read_csv(sys.argv[1], sep='\t', header=0, index_col=0)
df['Sex']=df['Sex'].apply(lambda x: 1 if x=="XY" else 2)
df.T.to_csv(sys.argv[2], sep='\t', header=True, index=True)
