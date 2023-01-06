import pandas as pd
import numpy as np

chunk = pd.read_csv("all_lines_test.csv", chunksize=10000, dtype=str)

df = pd.concat(chunk)

sample = df.sample(10)

print(sample)

sample.to_csv("sample.csv",index=False)
sample.save()