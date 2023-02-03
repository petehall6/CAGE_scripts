import pandas as pd


csv_file = input("Enter csv name: ")+".csv"

#open large csv via chunking
print("Loading CSV.\n")
chunk = pd.read_csv(csv_file, chunksize=10000)
df = pd.concat(chunk)

print(df.head(3))
print(df.count())


