import pandas as pd
import numpy as np
import datetime

data_store = "copy_number_processed.h5"

df = pd.read_hdf(data_store)

print(df.head())





