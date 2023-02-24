import os
import pandas as pd
import numpy as np
import openpyxl as xl
from pathlib import PureWindowsPath
from glob import glob

target_excel = "new_xlsx_1.xlsx"
dir = os.getcwd()
converted_dir = PureWindowsPath(r"Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\ngs_amplicon_count\converted")
ngs_folder_dest = PureWindowsPath(r"Z:\ResearchHome\Groups\millergrp\home\common\NGS\NGS 2022")
#count column b until not empty

def get_maxRows():
#loop through all files
    os.chdir(converted_dir)
    max_rows=0
    file_num=0
    print(os.getcwd())
    for file in os.listdir(os.getcwd()):
        f = os.path.join(converted_dir,file)
        if f.endswith(".xlsx"):
            try:
                df = pd.read_excel(f, engine="openpyxl")
                df = df.iloc[2]
                max_rows = max_rows + df.count()
                file_num = file_num + 1
            except:
                print(f"File failed: {f}")

    print(f"Number of amplicons: {max_rows}")
    print(f"Number of files counted: {file_num}")


def generate_xlsx():
    
    curr_dest = os.path.dirname(os.path.abspath(__file__))
    #os.chdir(curr_dest)

    target_folder = []
    print("Finding folders")
    #get a list of all folders in NGS folder
    for dir in os.scandir(ngs_folder_dest):
        if dir.is_dir():
            target_folder.append(dir.name)
    print(f"Number of target folders: {len(target_folder)}")       

    #parse through list of target_folders
    target_xls = []
    excel_count=0

    print("Finding excel sheet in folder")
    for dir in target_folder:
        os.chdir(os.path.join(ngs_folder_dest, dir))
        for f in os.scandir(os.path.join(ngs_folder_dest, dir)):
            if target_excel in f.name:
                excel_count = excel_count + 1
                target_xls.append(f)
            #else:
                #empty_folder.append(dir)

    print(f"Excel sheet count: {excel_count}")
    #print(f"Empty folders:  {empty_folder}")

    print("Converting xls to xlsx for parsing")
    #convert xls to xlsx
    failed_files=[]
    i=1
    os.chdir(curr_dest+"\\converted")
    print(os.getcwd())
    for f in target_xls:
        try:
            df = pd.read_excel(f) 
            with pd.ExcelWriter(f"new_xlsx_{i}.xlsx",engine="openpyxl") as writer:#keeps from throwing error with xlsx
                df.to_excel(writer)
                i = i+1
                del df
        except: #catches .xlsx files, next time just copy if already xlsx
            failed_file = os.path.abspath(f)
            failed_files.append(failed_file)
    print(f"Num of failed files: {len(failed_files)}")
    print(failed_files) 

#generate_xlsx()
get_maxRows()