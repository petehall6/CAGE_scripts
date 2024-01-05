import os
import argparse
import subprocess
import re
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
    return parser


def get_pam(casType):
    
    if casType == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif casType == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = casType

    return pam


def create_offinder_template(inputFile, pam, species):
    
    if species.lower() == 'm':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/mm10'
    
    if species.lower() == 'h':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/hg38'
    
    pam_df = pd.DataFrame([genome,pam])
    
    input_df = pd.read_csv(inputFile,sep='\t')    
    
    print(pam_df)
    print(input_df.head())
    
    output_df = pd.concat([pam_df, input_df], ignore_index=True)
    
    
    
    print(output_df.head())
    
    
    
    return






args = get_parser().parse_args()

inputFile = args.input
casType = args.cas
customType = args.custom
species = args.s

if customType != None:
    casType = customType
    
    
pam = get_pam(casType)

create_offinder_template(inputFile,pam,species)

#create_offinder_template(inputFile, casType)

