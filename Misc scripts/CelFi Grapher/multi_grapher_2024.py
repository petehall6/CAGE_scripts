import shutil
import pathlib
import os
import pandas as pd
from all_indels_barplot import find_csv, get_indels 

NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"

inputCSV = 'input_csv.csv'


csv_df = pd.read_csv(inputCSV, delimiter='\t', dtype=str)
crispy_list = csv_df.values.tolist()


for index in crispy_list:
    crispy,date,graph_type = index
    
    #adds the zero in front of month for Jan-Sept
    if len(date) != 6:
        date = "0"+date
    #run all_indels for each crispy file
    ngs_date_dir = os.path.join(NGS_DIR,date, "joined").replace("\\","/")
    
    #catches date
    if os.path.isdir(ngs_date_dir):
        #catches crispy
        print(f"\nDate found.  Attempting to run all_indels for {crispy}.")
        try:
            #go to NGS date folder
            os.chdir(ngs_date_dir)
            crispy_name = str(ngs_date_dir +"/"+ crispy + ".py")
            all_indel_crispy = crispy_name.replace(".py",f"_all_indels {graph_type} BP.py")
            
            input(crispy_name)
            input(all_indel_crispy)
            
            
            #make copy of crispy and name it all_indels
            shutil.copy(crispy_name,all_indel_crispy)
            
            #TODO get sample plates from csv in crispy folders.........
            
        except:
            None
    else:
        print(f"ERROR: all_indels for {crispy} not ran. \nCould not find NGS run date {date}. Please double check the csv.")
