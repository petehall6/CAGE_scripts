import pandas as pd

"""
Converts large copyNumber csv to hdf5 format for faster read times.
PMH 2023-01
    
"""

#copyNumber.csv is already scrubbed and lineage information is removed.

h5store = "copy_number_processed.h5"


#open large csv via chunking
print("Converting CSV to HDF.  This process can take a few minutes.\n")
chunk = pd.read_csv('copyNumber.csv', chunksize=10000)
df = pd.concat(chunk)
#df.astype(str)
df = df.rename(columns={"cell_line_display_name": "cell_line"})
print(df.head(3))
print(df.tail(3))


#convert to h5 file for quicker parsing
df.to_hdf(h5store, key='df', mode='w')
print("H5 file created\n\n")
