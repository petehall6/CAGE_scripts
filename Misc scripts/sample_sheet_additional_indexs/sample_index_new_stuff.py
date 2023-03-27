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



#make a giant list of index2 so that itll match index1 and i will save my sanity
all_index_list=[]

for index2 in index2_format:
    for index1 in index1_list:
        all_index_list.append(f'{index1}-{index2}')
#print(len(all_index_list))



#create sample id list
for plate in plates:#plates 80
    for well in well_num:#wells 96
        sample_ids.append(f"Miller-{plate}-{well}")
                
#print("Completed")
#print(len(sample_ids))
#print(sample_ids[-1])
    

final=[]

try:
    for i in range(len(sample_ids)):
        final.append(f"{sample_ids[i]}-{all_index_list[i]}")
    print(f"Length of final {len(final)}")
    print("Completed")
except:
    print(len(final))
    print(final[-1])
    print("failed")

df = pd.DataFrame(final)
df.to_csv('final_out.csv', index=False, header=False)