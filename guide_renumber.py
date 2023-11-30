import openpyxl as pxl
import numpy as np
import pandas as pd

#Be sure to sort the genes into by alphabetical inside of the csv
#This program assumes that the gene only appears once in the list
# Format the csv as follows

#Seq   #Gene  #Guide


#TODO Change the input and output names
input_csv ='lib109.csv' 
output_excel = 'lib109_guide_names.xlsx'


input_df = pd.read_csv(input_csv)

gene_list = input_df['Gene'].values.tolist()

print(input_df.head(5))
#print(gene_list)
guide_name_list = []
cur_gene =""
prev_gene = 'Dysf'
g_num = 1

for gene in gene_list:
    cur_gene = gene
    if cur_gene == prev_gene:
        guide_name_list.append(gene + '.g' + str(g_num))
        prev_gene = gene
        g_num +=1
    else:
        g_num = 1
        guide_name_list.append(gene + '.g' + str(g_num))
        prev_gene = gene
        g_num +=1

guide_name_df = input_df.assign(Guide=guide_name_list)

print(guide_name_df.head(5))

guide_name_df.to_excel(output_excel)

print('Complete')