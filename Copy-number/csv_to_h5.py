import pandas as pd
import numpy as np
import datetime

#copyNumber.csv is already scrubbed and lineage information is removed.

h5store = "copy_number_processed.h5"


#open large csv via chunking
chunk = pd.read_csv('copyNumber.csv', chunksize=10000)
df = pd.concat(chunk)
df.astype(str)
#print(df.tail())
df = df.rename(columns={"cell_line_display_name": "cell_line"})
print(df.head(3))


#convert to h5 file for quicker parsing
df.to_hdf(h5store,key='df', mode='w')
print("H5 file created\n\n")
