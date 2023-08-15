import pandas as pd
from Bio.Seq import Seq


fwd_backbone = 'CACCG'
rev_backbone = 'AAAC'
rev_suffix = 'C'


#get input from excel sheet
wb = pd.read_excel('input.xlsx')

#set starting index to 1
wb.index += 1

#setup temporary df for renaming guides and getting reverse compliment
out_columns = ['Forward Name', 'Forward Sequence', 'Reverse Name', 'Reverse Sequence']
tmp = pd.DataFrame(columns=out_columns)
tmp  = pd.concat([wb, tmp])

output_fwd_names = [name+'.F' for name in tmp['Guide Name']]
output_rev_names = [name+'.R' for name in tmp['Guide Name']]

fwd_seq = [Seq(seq) for seq in tmp['Seqeunce']]
rev_seq = [seq.reverse_complement() for seq in fwd_seq]

fwd_seq_oligo = [str(fwd_backbone + seq) for seq in fwd_seq]
rev_seq_oligo = [str(rev_backbone + seq + rev_suffix) for seq in rev_seq]

tmp['Forward Name'] = output_fwd_names
tmp['Reverse Name'] = output_rev_names
tmp['Forward Sequence'] = fwd_seq_oligo
tmp['Reverse Sequence'] = rev_seq_oligo

#quick way to merge 4 columns into two
grouped_name = []
grouped_seq = []
for index, row in tmp.iterrows():
    grouped_name.append(row['Forward Name'])
    grouped_name.append(row['Reverse Name'])
    grouped_seq.append(row['Forward Sequence'])
    grouped_seq.append(row['Reverse Sequence'])


out_wb = pd.DataFrame(columns=['Name','Sequence'])
out_wb['Name'] = grouped_name
out_wb['Sequence'] = grouped_seq

print(out_wb)

out_wb.to_excel('output.xlsx', index=False,index_label=False)