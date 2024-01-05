import os
import argparse
import subprocess
import re
import time
import pandas as pd

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input',
        type=str,
        help='File with raw sequencing reads')
    parser.add_argument('-s',
        type=str,
        help='Species')
    parser.add_argument('-cas',
        type=str,
        required=False,
        help='cas')
    parser.add_argument('-custom',
        type=str,
        required=False,
        help='Manually enter custom PAM')
    parser.add_argument('-n',
        type=int,
        required=False,
        help='Number of guides per gene. Defaults to 4 if nothing specified')
    
    
    return parser

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
    print(output_df.head())

def offinder():

    input_guides = 'columns_combined_for_offinder.txt'

    #have to load primer3 to get it to work.
    os.system('module load primer3')
    os.system('module load gcc')
    os.system('module load cas-offinder')

    offinder_call = f'bsub -P Cas_offinder -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 8000 -o out_CasOffinder.out -e err_CasOffinder.err -J cas_offinder "hostname & cas-offinder {input_guides} G CasOffinder_Results.txt"'
                    

    Long_sout = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')    
    jid_Long9 = re.split("<|>", Long_sout)[1]
    if jid_Long9.isnumeric():
        print(f"JobID-{jid_Long9} for cas_offinder has been successfully submitted. Please wait.")
        print("System will notify when the job is complete")
    else:
        print("Failed in submitting job for cas-offinder OTA for Cas9Long")
    os.system(f"bwait -w 'ended(cas_offinder)'")
    print("Offinder complete")
    print("Taking top picks")

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

#*passing arguments block
'''
args = get_parser().parse_args()

inputFile = args.input
casType = args.cas
customType = args.custom
species = args.s
grna_num = args.n

'''

#*debug block
inputFile = 'test_offinder_input.txt'
casType = '9'
customType = None
species = 'h'
grna_num = 5 



if customType != None:
    casType = customType
if grna_num == None:
    grna_num = 4


pam = get_pam(casType)

create_offinder_template(inputFile,pam,species)
offinder()
#clean_files()


