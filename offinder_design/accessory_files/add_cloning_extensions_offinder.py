import pandas as pd
import numpy as np
import os
from pathlib import Path


def add_extension(lib_name, desired_primer, input_filename):

    LIB_NAME = lib_name # CHANGE ME

    DESIRED_PRIMER = desired_primer # CHANGE ME (IF NECESSARY) PC-155 for gibson
    INPUT_FILENAME = input_filename # CHANGE ME (IF NECESSARY)

    INPUT_DIR = Path.cwd() 

    OUTPUT_FILENAME = f'{LIB_NAME}-{DESIRED_PRIMER}_twist.tsv'

    TWIST_PREFIX = f'{LIB_NAME}_'

    LIB_AND_PRIM_DIR = os.path.abspath(r'/research_jude/rgs01_jude/groups/millergrp/home/common/Screens')
    LIB_AND_PRIM_XLSX = os.path.join(LIB_AND_PRIM_DIR, 'lib_vectors_and_primers.xlsx')

    ext_df = pd.read_excel(LIB_AND_PRIM_XLSX, sheet_name='extension_for_TWIST', engine='openpyxl') # Extensions DF
    # mrp_df = pd.read_excel(LIB_AND_PRIM_XLSX, sheet_name='cloning_mrp', index_col=0) # "Most recent Primers" df
    primers_df = pd.read_excel(LIB_AND_PRIM_XLSX, sheet_name='extension_for_TWIST', index_col='primer_set_name', engine='openpyxl')
    primer_info = primers_df[primers_df.index == DESIRED_PRIMER].head(1).squeeze()
    primer_info = primer_info.fillna('')

    # Cloning oligo sequence is:
    #   5' primer + 5' BsmbI + gRNA (no PAM) + 3' BsmbI + 3' primer
    five_primer = primer_info["5_extension"]
    five_bsmbi = primer_info["5_BsmBI"]
    three_bsmbi = primer_info["3_BsmBI"]
    three_primer = primer_info["3_extension"]

    input_df = pd.read_excel(INPUT_FILENAME, engine='openpyxl')


    input_df['no_pam'] = input_df['full_gRNA'].apply(lambda seq: seq[:20]) # Make a seq w/o PAM to join Broad info
    input_df['oligo_w_extensions'] = five_primer + five_bsmbi + input_df['no_pam'] + three_bsmbi + three_primer

    print(f'OLIGO TEMPLATE: {five_primer} {five_bsmbi} [ 20bp gRNA seq ] {three_bsmbi} {three_primer}')

    oligo_qc = {
        'avg': np.average(input_df['oligo_w_extensions'].map(len)),
        'stDev': np.std(input_df['oligo_w_extensions'].map(len))
    }
    no_pam_qc = {
        'avg': np.average(input_df['no_pam'].map(len)),
        'stDev': np.std(input_df['no_pam'].map(len))
    }

    print(f"""
    oligo QC:
        Avg   - {oligo_qc['avg']}
        StDev - {oligo_qc['stDev']}

    no_pam_QC:
        Avg   - {no_pam_qc['avg']}
        StDev - {no_pam_qc['stDev']}
    """)


    input_df['twist_name'] = TWIST_PREFIX + input_df['Name']

    output_df = input_df[['twist_name', 'oligo_w_extensions']]
    output_df.columns = ['name', 'seq']

    output_df.to_csv(OUTPUT_FILENAME, header=False, index=False, sep='\t')

    print('Done!')
    print('Output File:')
    print(f'    {OUTPUT_FILENAME}')
