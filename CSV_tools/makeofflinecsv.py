##convert csv from wide to long.  Sqlite3 column limit = 2000.
import pandas as pd


data_store = "copyNumSql.csv"
long_csv = "copy_number_long.csv"
copy_db = "copy_num.db"

#open large csv via chunking
print("Chunking CSV.  Could take a while\n")
chunk = pd.read_csv(data_store, chunksize=1000)
df = pd.concat(chunk)
df = df.rename(columns={"cell_line_display_name": "cell_name"})
#will convert to int
df = df.round()
print("Chunk complete\n")

df2 = df.sample(100)

df2.to_csv('offlinecsv-wide.csv')