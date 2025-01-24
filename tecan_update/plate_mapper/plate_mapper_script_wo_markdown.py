import os
import pandas as pd
import math

def mapped_reactions(reactions):
    
    if reactions % 2 == 0:
        return 0
    else:
        return 1



# Load the Excel file 'input.xlsx' into a DataFrame
input_df = pd.read_excel('input.xlsx', engine='openpyxl')

input_df.index =input_df.index + 1

#Calculate amount of DNA in sample
input_df['DNA ug'] = input_df['Volume'] * (input_df['Concentration nanograms'] / 1000)

#Calculate coverage
input_df['ug of DNA needed'] = (input_df['number of guides'].astype(float) * 0.006 * input_df['Coverage']) / 1000

#Round up to nearest whole number
input_df['number of reactions'] = (((input_df['Volume'] * input_df['Concentration nanograms']) / 1000 ) / input_df['ug of DNA needed']).apply(lambda x: math.ceil(x))

#Add addtional reaction to odd numbers so that plate maps are evenly spaced
input_df['empty reactions'] = input_df['number of reactions'].apply(mapped_reactions)

#TODO
#Will need to limit each reaction to 10ug/rxn then find minimum number of reactions to reach desired coverage

#has the extra reaction
rxn_df = input_df[['Sample','number of reactions', 'empty reactions']].copy()

rxn_df

rxn_list = rxn_df.values.tolist()

print(rxn_list)

sample_rxn_list = []
for sample in rxn_list:
    for rxn in range(sample[1]):
        sample_rxn_list.append(sample[0])
    if sample[2] == 1:
        sample_rxn_list.append('empty')

print(sample_rxn_list)

well_cols = list(range(1,13))
well_rows = ['A','B','C','D','E','F','G','H']

well_dict = {}

for row in well_rows:
    for col in well_cols:
        well_dict[row + str(col)] ='empty'


for i, well in enumerate(well_dict):
    try:
        well_dict[well] = sample_rxn_list[i]
    except:
        well_dict[well] = 'empty'

print(well_dict)



#generating a dict of rows ord() returns the unicode value of a character, chr() returns the character of a unicode value
rows_dict = {}
for i in range(ord('A'), ord('I')):
    row = chr(i)
    rows_dict[row] = []
    
#filling the rows with the well_dict values
for well in well_dict.items():
    for row in rows_dict.items():
        if well[0][0] == row[0]:
            row[1].append(well[1])
            
#converting the dict to a dataframe
plate_layout = pd.DataFrame(rows_dict)
#transposing the dataframe, going from long to wide
plate_layout = plate_layout.T

#renumber columns to match plate
plate_layout.columns = list(range(1,13))

print(plate_layout)

plate_layout.to_excel('plate_layout.xlsx', engine='openpyxl')


