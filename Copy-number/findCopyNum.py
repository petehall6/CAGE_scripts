import pandas as pd
import numpy as np
import datetime

data_store = "copy_number_processed.h5"

df = pd.read_hdf(data_store)

#strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
#stripping here does not impact performance
df.astype(str)
df.columns= df.columns.str.replace("-","")

cell_line = "143B"
gene = "ADA"

print(df.loc[cell_line,gene])
print(df.get(gene, default='Cell not found'))





