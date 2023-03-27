from openpyxl import Workbook
import openpyxl
import os
import glob
import shutil
import pandas as pd



def convertoExcel():
    cur_dir = os.getcwd()
    new_excel_list=[]
    mmx_full = glob.glob(os.path.join(cur_dir,'*mmx_full.tsv'))
    mmx_gen = glob.glob(os.path.join(cur_dir,'*mmx_gen.tsv'))
    mmx_simp= glob.glob(os.path.join(cur_dir,'*mmx_simplified.tsv'))
    mmx_input =glob.glob(os.path.join(cur_dir,'input.tsv'))

    mmx_full=str(mmx_full[0])
    mmx_gen = str(mmx_gen[0])
    mmx_simp = str(mmx_simp[0])
    mmx_input = str(mmx_input[0])

    input_tsv_list=[mmx_gen,mmx_full,mmx_simp,mmx_input]
    print(input_tsv_list)
    #access TSV cellls
    for tsv in input_tsv_list:
        name = tsv[-8:].strip("_").replace(".tsv",".xlsx")
        if name =='fied.tsv':
            name ='simp.tsv'
        print(name)
        df = pd.read_excel(tsv, sep="\t")
        with pd.ExcelWriter(name,engine='openpyxl') as writer:
            df.to_excel(writer)
            new_excel_list.append('mmx_full.xlsx')
            writer.close()
    print(new_excel_list)
    #if folder not exists
    os.mkdir("new_excel")
    os.chdir(os.path.join(cur_dir,"new_excel"))

    new_dir = os.getcwd()
    for file in new_excel_list:
        shutil.move(os.path.join(cur_dir,file),new_dir)

#have TSV need to change to .xslx


convertoExcel() 




