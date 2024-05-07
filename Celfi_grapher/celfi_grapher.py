import pandas as pd
import matplotlib as plt
import os
import shutil
from pathlib import Path
import glob
import math


inputCSV = "input.csv"
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
csv_dir = os.path.join(os.getcwd(),'csv files')
print(csv_dir)

def find_csv():
    #get input from csv
    input_df = pd.read_csv(inputCSV)

    #puts columns into list to containerize the data eaiser
    input_projects = input_df.values.tolist()
    print("Scanning for csv")
    
    for target_proj in input_projects:
        
        cage_proj, ngs_date = target_proj
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined')
        os.chdir(joined_dir)
        
        #search joined folder for all_indel
        for file in glob.glob(f"{cage_proj}**\*all_indels.csv",recursive=True):
            if cage_proj in file:
                shutil.copy(file, csv_dir)
    
    print("CSV files copied over.")

def make_plot():
    
    os.chdir(csv_dir)
    
    columns = ['Sample','Out-of-frame']
    
    for file in os.scandir(os.getcwd()):
        oof_df = pd.read_csv(file,usecols=columns)
        oof_df.index = oof_df.index + 1
        print(oof_df.head(10))
        
        comparison_number = int(input("How Many CelFi Scores: "))
        
        comparison_counter = 0
        comparisons_list = []
        while comparison_counter < comparison_number*2: #doubled since each score takes 2 data points
            curr_comparison = math.ceil((comparison_counter+1)/2)
            
            comparision_choice = int(input(f"Enter ROW number for comparison {math.ceil((comparison_counter+1)/2)} : "))
            
            #do initial and final time points
            #dont bother with trying to be smart and having the loop handle all the guessing
            #will need a loop inside this one to handle inidividual comparison then have a counter to check when all comparisons are finished
            
            comparison_counter +=1
        
        print(f"Here are the comparions: {comparisons_list}")
        
    
    return None 

#find_csv()

make_plot()