import csv
import pandas as pd

index1_list = []



sample_ids=[]



well_num=[]
sample='Miller'
plates=[]
index2_raw ='CAGCGTGATCACACCAGTTGCACGACATTAGTGTAGCTAGTCTGTGCATCAGGACGGTTATTAACTATGAACCCTAAGAATCCGGGCTGCTACCTTTCTTAAGTCAGGATACTGTATGTC'
#called once to get rid of giant 385 line string. Copied indexes from Excel and formatted them. Then saved as external csv to work with later

""" #create index list
index_format=index.split('\n')
#dump empty newline
index_format=index_format[0:384]
print(len(index_format))
print(index_format)
#save csv of formatted index to clean up script.  Will call import csv list later
df = pd.DataFrame(index_format)
df.to_csv('index1_out.csv', index=False, header=False) """

index1_list=[]
index1_file=open("index1_out.csv", "r")
index1_format=list(csv.reader(index1_file, delimiter=","))
index1_file.close

for index in index1_format:
    #converts from a giant list of list to single indexs
    a= ''.join(index)
    index1_list.append(a)
    
#print(index1_list)





#create well_number
for letter in range(65,73):
    row=chr(letter)
    for num in range(1,13):
       column=str(num)
       well_loc=str(row)+str(num)
       well_num.append(well_loc)
#print(well_num) 

#create plate#
for sample in range(81,161):
    plates.append(f'Plate{sample}')   
#print(plates)
    

#create index2 list
index2_format=[]
n=6
for i in range(0,len(index2_raw),n):
    index2_format.append(index2_raw[i:i+n])
#print(index2_format)
#print(len(index2_format))

giant_index2=[]





print(len(index2_format)) 
print(len(giant_index2))   



""" 
print(plates)
print(f"Well List len {len(well_num)}")
print(f"Plate list len {len(plates)}")
print(f"Index1 List len {len(index1_list)}")
print(f"Index2 List len {len(index2_format)}")
"""
"""

index_set=0

try:
    for i in range(80):#plates 80
        index_set=0 
        for ii in range(20):#index 2 20
            index_set = index_set+1
            for iiii in range(96):#wells 96
                if index_set == 1:
                    sample_ids.append(f"Miller-{plates[i]}-{well_num[iiii]}-{index1_list[iiii]}-{index2_format[ii]}")
                if index_set == 2:
                    sample_ids.append(f"Miller-{plates[i]}-{well_num[iiii]}-{index1_list[iiii+96]}-{index2_format[ii]}")
                if index_set == 3:
                    sample_ids.append(f"Miller-{plates[i]}-{well_num[iiii]}-{index1_list[iiii+192]}-{index2_format[ii]}")
                if index_set == 4:
                    sample_ids.append(f"Miller-{plates[i]}-{well_num[iiii]}-{index1_list[iiii+288]}-{index2_format[ii]}")
                    
    print("Completed")
    print(len(sample_ids))
    print(sample_ids[-1])
except IndexError:
    print(len(sample_ids))        
    print(sample_ids)
    print(IndexError) """

