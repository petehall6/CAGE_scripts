import pandas as pd
import numpy as np
import datetime

#copyNumber.csv is already scrubbed and lineage information is removed.

h5store = "copy_number_processed.h5"


#open large csv via chunking
chunk = pd.read_csv('copyNumber.csv', chunksize=10000)
df = pd.concat(chunk)
df = df.rename(columns={"cell_line_display_name": "cell_line"})
df.set_index("cell_line", drop=False)
print(df.tail())
#convert to h5 file for quicker parsing
df.to_hdf(h5store,key='df', mode='w')
print("H5 file created\n\n")

df = pd.read_hdf(h5store,"df")

print(df.tail())








