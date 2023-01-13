import pandas as pd
import numpy as np
from dask import dataframe as ddf
import sqlite3
import psycopg2


#copyNumber.csv is already scrubbed and lineage information is removed.


#open large csv via chunking
chunk = pd.read_csv('copyNumber.csv', chunksize=1000)

df = pd.concat(chunk)

#Check for cell information 

print(df.head())

myCell = df.iloc['143B']
print(myCell)


