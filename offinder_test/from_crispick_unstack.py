import os
import argparse
import subprocess
import re
import time
import pandas as pd
import numpy as np

inputFile = "lib126grna-designs.txt"
casType = '9'
customType = None
species = 'm'


#*Need to add in Neg controls back into library
#*Need to sort guides alphabetically then assign guide names before sending to offinder...or after taking top hits and then adding g#s.


def get_pam(casType):
    
    if casType == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif casType == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = casType

    return pam

def combine_columns(rows):
    
    if pam == 'NNNNNNNNNNNNNNNNNNNNNGG':
        if str(rows['Target Gene Symbol']) != '':
            return str(rows['sgRNA Sequence']) + str(rows['PAM Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])
        
    elif pam =='TTTVNNNNNNNNNNNNNNNNNNNNN':
        if str(rows['Target Gene Symbol']) != '':
            return str(rows['PAM Sequence']) + str(rows['sgRNA Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])

def format_crispick(inputFile, pam, species):
    
    if species.lower() == 'm':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/mm10'
    
    if species.lower() == 'h':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/hg38'
        
    pam_df = pd.DataFrame([[genome],[pam]],columns=['Combined'])

    input_df = pd.read_csv(inputFile, sep='\t', usecols=['Target Gene Symbol','sgRNA Sequence','PAM Sequence'])
    
    input_df['Combined'] = input_df.apply(combine_columns, axis=1)

    input_df = input_df.drop(columns=['Target Gene Symbol','sgRNA Sequence'])
    input_df.drop(columns=['PAM Sequence'], inplace=True)

    output_df = pd.concat([pam_df,input_df],ignore_index=True)

    output_df.to_csv("columns_combined_for_offinder.txt",sep='\t',header=False,index=False)    
        
    #make seperate df for NTCs

def offinder():

    input_guides = 'columns_combined_for_offinder.txt'

    offinder_call = f'bsub -P Cas_offinder -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 8000 -o run_summary.txt -e error_summary.txt -J cas_offinder "hostname & cas-offinder {input_guides} G CasOffinder_Results.txt"'
                    
    Long_sout = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')    
    jid_Long9 = re.split("<|>", Long_sout)[1]
    if jid_Long9.isnumeric():
        print(f"JobID-{jid_Long9} for cas_offinder has been successfully submitted. Please wait.")
        print("System will notify when the job is complete")
    else:
        print("Failed in submitting job for cas-offinder OTA for Cas9Long")
    os.system(f"bwait -w 'ended(cas_offinder)'")
    print("Jobs done")

def create_ota_table():
    input_file = 'CasOffinder_Results.txt'

    #summary frame serves as the template for the completed mismatch dataframe
    #Count_alignment for Cas 9Long
    tmp_df =pd.read_csv(input_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None, index_col=False)

    tmp_df.drop(['chrom', "strand", 'result'],axis=1,inplace=True)
    
    count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()

    count_df.columns = count_df.columns.droplevel(1)
    count_df.fillna(0, inplace=True)
    
    count_df.columns = ['Long_0','Long_1','Long_2','Long_3']
    
    count_df['Long_1']=count_df['Long_1']+count_df['Long_0']
    count_df['Long_2']=count_df['Long_2']+count_df['Long_1']+count_df['Long_0']
    count_df['Long_3']=count_df['Long_3']+count_df['Long_2']+count_df['Long_1']+count_df['Long_0']
    
    count_df.reset_index(inplace=True)
    count_df.insert(1,"Gene",np.nan)
    
    
    
    #count_df['gRNA seq'] = count_df.index
    #count_df.index
    
    print(tmp_df['name'])
    print(count_df.shape)
    print("*"*20)
    print(count_df.head(10))
    print("*"*50)

    count_df.to_csv("OTA_scores_sorted.csv")

def clean_files():
    #waits to make sure all files are done writing before going to next step
    time.sleep(3)    
    try:
        os.remove("out_CasOffinder.out")
        os.remove("err_CasOffinder.err")
        os.remove("columns_combined_for_offinder.txt")
        print("files cleaned")
    except:
        print("didnt delete anything.")

if customType != None:
    casType = customType

pam = get_pam(casType)

#format_crispick(inputFile,pam,species)

#offinder()

create_ota_table()


#clean_files()