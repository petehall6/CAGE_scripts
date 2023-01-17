import pandas as pd
import numpy as np
import time

#read in HDF5 file from csv_to_h5 script
data_store = "copy_number_processed.h5"

#for print outline of copy number
top_btm_boarder = "*"*33
line_boarder = "*"+" "*31 + "*"

#read_time_start = time.time()
df = pd.read_hdf(data_store)
#read_time_end = time.time()
#print("File Time: ", read_time_end-read_time_start)

#strips hyphens for easier parsing.  Trying to remove hyphen in the hdf5 will give an unique index error
#stripping here does not impact performance

#print(df.head(3))

#TODO is this even needed? cast as strings

#df.apply(str)

#strip_start_time = time.time()
df.columns= df.columns.str.replace("-","")
#strip_end_time = time.time()
#print("Strip time: ", strip_end_time-strip_start_time)


cell_line = input("Enter cell line: ").upper()
gene = input("Enter gene: ").upper()

#wont work because index is not cell_line....what happens if index is set to cell_line
print(df.loc[df[cell_line], gene])





















#get cell line and gene information for query



try:
    copyNumber = str(df.loc[cell_line,gene])
    
    print(top_btm_boarder)
    print(line_boarder)
    print("\tNumber of copies: \t",copyNumber[0])
    print(line_boarder)
    print(top_btm_boarder)
    
except KeyError:
    print("Problem finding cell line: {} or gene: {}".format(cell_line, gene))

#possible matches
#Gene matches Work



#Cell lines not working

#filters take a function so lambda was used to create anyonmous functions
cell_matches = list(filter(lambda x: cell_line in x, df['cell_line']))
gene_matches = list(filter(lambda x: gene in x, df.columns))

























"""
for (i,j) in enumerate(cell_matches, start=1):
    print(i,j)

cell_pick = int(input("Enter cell choice number: "))



for (i,j) in enumerate(gene_matches, start=1):
    print(i,j)

gene_pick = int(input("Enter gene choice number: "))

cell_choice = cell_matches[cell_pick-1]
gene_choice = gene_matches[gene_pick-1]

print("Cell choice: ",cell_matches[cell_pick-1])
print("Gene choice: ",gene_matches[gene_pick-1])

confirm = input("Press 'y' to confirm: ").upper()

if confirm == 'Y':
    copyNumber = str(df.loc[cell_choice,gene_choice])
    print(copyNumber) """





#print(df.get(gene, default='Cell not found'))





