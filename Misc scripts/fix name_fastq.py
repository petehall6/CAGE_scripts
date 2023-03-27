import os
import fnmatch

file_list=[]

data_dir = os.getcwd()

a1_indexs = ["GGTCAGAGCA","ATAGAGCCAG","GCTGAGGGAA","TCTTATTTAG"]

for file in os.listdir(data_dir):
    if file.endswith(".fastq") and "B9" in file:
        file_list.append(file)

print(len(file_list))

""" for name in file_list:
    new_name = name.replace('A10', 'A01')
    os.rename(name,new_name)
     """
    
for file in os.listdir(data_dir):
    if file.endswith(".fastq") and "B9" in file:
        new_name = file.replace('B9', 'B09')
        os.rename(file,new_name)
    if file.endswith(".fastq") and "B8" in file:
        new_name = file.replace('B8', 'B08')
        os.rename(file,new_name)


