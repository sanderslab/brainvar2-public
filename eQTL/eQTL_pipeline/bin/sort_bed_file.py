#!/usr/bin/env python3
import sys
import pandas as pd

bed_file = sys.argv[1]
sorted_bed_file = sys.argv[2]

bed = pd.read_csv(bed_file, sep='\t')
x=bed.iloc[:, 4:]
first_four = bed.iloc[:, :4]
sorted=pd.concat((first_four, x.sort_index(axis=1)), axis=1)
sorted.to_csv(sorted_bed_file, sep ='\t', index=False)
