import pandas as pd
import os

#TODO why doesnt this run on the cluster?


#* Change these as needed.  Don't touch anything else.
# run on local machine


input_file = 'CasOffinder_Results.txt'
output_file = 'custom_PAM_for_PC2.csv'
long_mismatch_num = 3



os.chdir("input")
#summary frame serves as the template for the completed mismatch dataframe
#Count_alignment for Cas 9Long
df =pd.read_csv(input_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None)
df.sort_values(by=['gRNA', 'mismatch'])

summary_frame = pd.DataFrame(columns=["name","gRNA","Long_0","Long_1","Long_2","Long_3"])
summary_frame["name"] = df['name']
summary_frame["gRNA"] = df['gRNA']
summary_frame.fillna(0,inplace=True)


for iMismatch in range(0, long_mismatch_num+1):
    for igRNA in set(df['gRNA']):       
        tmp=df.loc[(df['gRNA'] == igRNA) & (df['mismatch'] <= iMismatch)]
        seq_name=summary_frame[summary_frame['gRNA']==igRNA].index.tolist()
        column='Long_'+str(iMismatch)
        summary_frame.loc[seq_name, column] = len(tmp)
os.chdir("../output")
summary_frame.to_csv(output_file) 
