import os
import subprocess
import pandas as pd
import time
import re
import numpy as np

#*Update Values as needed
#*output should always end in .txt

inputFile = 'input_template.txt'
outputFile = 'gene_list_result_test.txt' #will output txt files
casType = '9' #or 12
species = 'h' #'m' for mouse, 'h' for human.  Currently only supported species but can adjust as needed


'''
Make sure your input file is saved as a .txt and formatted as follows.  Note it includes the pam ie NGG:
Sequence    Name
GGGAGGCGGTGGCATCAGTTCAGANGG    .g1

The columns dont have to align as long as they seperated by a tab.  
It's easiest to create the .txt in Excel and make Column A Sequence and Column B Name.
There should be a demo file named test_offinder_input.txt in this folder. 
'''

def get_pam(casType):
    
    if casType == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif casType == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = casType

    return pam

def combine_columns(rows):
    
    return rows['Sequence']+ ' 3 ' +rows['Name']

def create_offinder_template(inputFile, pam, species):
    
    if species.lower() == 'm':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/mm10'
    
    if species.lower() == 'h':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/hg38'
        
    pam_df = pd.DataFrame([[genome],[pam]],columns=['Combined'])
    
    input_df = pd.read_csv(inputFile,sep='\t')    
    #concats the sequence, mismatch and name all on one line
    input_df['Combined'] = input_df.apply(combine_columns, axis=1)
    input_df = input_df.drop(columns=['Sequence','Name'])


    output_df = pd.concat([pam_df,input_df],ignore_index=True)

    output_df.to_csv("columns_combined_for_offinder.txt",sep='\t',header=False,index=False)
    
    print("This is going into Offinder")    
    print("\n\n")
    print(output_df)

def offinder():

    #have to load primer3 to get it to work.
    offinder_call = 'bsub -P Cas_offinder -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 8000 -o offinder_job_summary.txt -e error_log.txt -J Cas_offinder cas-offinder "hostname & cas_offinder columns_combined_for_offinder.txt G offinder_results.txt"'
    
    Long_sout = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')
    jid_Long9 = re.split("<|>", Long_sout)[1]
    if jid_Long9.isnumeric():
        print(f"JobID-{jid_Long9} for cas_offinder has been successfully submitted. Please wait.")
        print("System will notify when the job is complete")
    else:
        print("Failed in submitting job for cas-offinder OTA for Cas9Long")
    os.system(f"bwait -w 'ended(Cas_offinder)'")
    print("Offinder complete")

def ota_scores():
    input_file = 'offinder_results.txt'

    #removes decimal place
    pd.set_option('precision',0)
    
    #summary frame serves as the template for the completed mismatch dataframe
    #Count_alignment for Cas 9Long
    tmp_df =pd.read_csv(input_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None)
    tmp_df.drop(['chrom', "strand","name", 'result'],axis=1,inplace=True)
    count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()
    count_df.columns = count_df.columns.droplevel(0)
    count_df.reset_index(0)
    count_df.fillna(0, inplace=True)
    
    #checks for missing columns and fills in the blank
    columns_needed = set([0,1,2,3])
    column_set = set(count_df.columns)
    missing_cols = list(sorted(columns_needed - column_set))
    for column in missing_cols:
        count_df.insert(loc=column, column=column, value=[0 for i in range(count_df.shape[0])])
    
    count_df.columns = ['Long_0','Long_1','Long_2','Long_3'] 
    #add Long_3 first so as not to change the value of Long_0,Long_1 and Long_2 before adding up to Long_3
    count_df['Long_3'] = count_df.iloc[:,0:4].sum(axis=1)
    count_df['Long_2'] = count_df.iloc[:,0:3].sum(axis=1)
    count_df['Long_1'] = count_df.iloc[:,0:2].sum(axis=1)

    count_df.sort_values(by=["Long_0", "Long_1", "Long_2","Long_3","gRNA"],inplace=True, ascending=True)
    
    print(count_df)
    
    #save to file
    count_df.to_csv(outputFile,header=True,sep='\t')

    print("******** OTA scores compiled ***************")



def create_ota_table():
    casoff_file = 'offinder_results.txt'

    #summary frame serves as the template for the completed mismatch dataframe
    tmp_df =pd.read_csv(casoff_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None, index_col=False)

    tmp_df.drop(['chrom', "strand", 'result'],axis=1,inplace=True)
    
    #get gene names now to merge later
    gene_df = tmp_df[['gRNA','name']]
    
    #groups df by the position, mismatch and gRNA seq.  Counts the number of mismatches per position and creates
    #multi-index df that lists the positions as columns for counting across the series        
    count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()

    #dropping one of the indexes for a flatter df
    count_df.columns = count_df.columns.droplevel(1)
    count_df.fillna(0, inplace=True)
    
    count_df.columns = ['Long_0','Long_1','Long_2','Long_3']
    
    #counts across the series of positional mismatches and renames them Long_0 etc
    count_df['Long_3'] = count_df.iloc[:,0:4].sum(axis=1)
    count_df['Long_2'] = count_df.iloc[:,0:3].sum(axis=1)
    count_df['Long_1'] = count_df.iloc[:,0:2].sum(axis=1)
    
    #reset index to get rid of any multi-index shannigans from above step and merge with gene names
    #gets rid of duplicatese.  count_df will not have any NTC's.  
    count_df.reset_index(inplace=True)
    count_df = count_df.merge(gene_df, on="gRNA", how="right",right_index=False).drop_duplicates()
    #ignore index to reset it
    sorted_df = count_df[['gRNA','name','Long_0','Long_1','Long_2','Long_3']].sort_values(by=['name','Long_0','Long_1','Long_2','Long_3'],ignore_index=True)

    #change capatilization of gene name before addding NTC to avoid switchcasing NTC
    if species == 'h':
        sorted_df['name'] = sorted_df['name'].str.upper()
    if species =='m':
        sorted_df['name'] = sorted_df['name'].str.capitalize()
    
    #add in NTC's from CRISPICK file, adds NTCs to bototm of list 
    #ntc_df = get_non_target_control()
    #sorted_df.sort_values(by=['Long_0'],ascending=True,inplace=True)    
    #sorted_df = pd.concat([sorted_df,ntc_df])
    
    #rearrange and rename columns to fit downstream (library draft) workflow
    arranged_columns = ['name','gRNA','Long_0','Long_1','Long_2','Long_3']
    sorted_df = sorted_df[arranged_columns]
    sorted_df = sorted_df.rename(columns={'name':'Name','gRNA':'full_gRNA'})

    #sorted_df = (sorted_df.pipe(add_guide_numbers))
    
    print(sorted_df.head())
    
    #sorted_df.to_excel(sortedFileName,index=False)
    

def clean_files():
    #waits to make sure all files are done writing before going to next step
    time.sleep(3)    
    try:
        os.remove("error_log.txt")
        os.remove("columns_combined_for_offinder.txt")
        os.remove("offinder_job_summary.txt")
        print("files cleaned")
    except:
        print("didnt delete anything.")

pam = get_pam(casType)

create_offinder_template(inputFile,pam,species)
offinder()
create_ota_table()
#ota_scores()
clean_files()