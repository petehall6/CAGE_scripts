import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import glob
import re
import numpy as np

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

def get_scores():
    def _pick_guides(comp):

            if comp ==1 :
                initial_time_pt_choice = 1
                final_time_pt_choice = 4
            else:
                initial_time_pt_choice = 5
                final_time_pt_choice = 8
        
        
            initial_time_pt_choice = int(input(f"initial time point for comparison {comp}: "))
            final_time_pt_choice = int(input(f"final time point for comparison {comp}: "))
            
            init_time_name = str(oof_df.loc[initial_time_pt_choice][0]).partition('g')
            init_time_guide = init_time_name[1]+init_time_name[2]
            init_time_oof = (oof_df.loc[initial_time_pt_choice][1])
            
            final_time_name = str(oof_df.loc[final_time_pt_choice][0]).partition('g')
            final_time_guide = final_time_name[1]+final_time_name[2]
            final_time_oof = (oof_df.loc[final_time_pt_choice][1])

            #checks to see if guides match.  if not it will restart
            if init_time_guide != final_time_guide:
                input("Guides do NOT match.  Press Enter to repick guides: ").lower()
                get_scores()
                
            guide_name = init_time_guide
            
            return guide_name, init_time_oof, final_time_oof
    
    fa_score_list=[]
    os.chdir(csv_dir)
    
    for csv in os.scandir(csv_dir):
        proj_name = csv.name.split("_")[0]
        print(proj_name)

        oof_columns = ['Sample','Out-of-frame']
        oof_df = pd.read_csv(csv,usecols=oof_columns)
        oof_df.index = oof_df.index+1

        print(oof_df)
        
        #comp_num=2
        
        comp_num = int(input("Number of Comparisons: "))
        
        
        
        for comp in range(1, comp_num+1):
            print(f"Comparions #: {comp}")
            guide_name, init_time_oof, final_time_oof = _pick_guides(comp)
            
            fa_score = round((final_time_oof/init_time_oof),2)
            print(f"\n\nFitness assay score for comparison {comp}: {fa_score}\n\n")
            fa_score_list.append([guide_name, fa_score])

        print(f"Fitness scores: {fa_score_list}")
        
    return fa_score_list

def graph_scores():
    fa_scores = get_scores()
    
    fa_score_df=pd.DataFrame(fa_scores, columns=["guide","fa_score"])
    
    y_max = fa_score_df['fa_score'].max()
    
    print(fa_score_df)
    
    fa_plot = fa_score_df.plot(kind="bar",
                               x="guide",
                               y='fa_score',
                               rot=0,legend=False,
                               title="Dependency Scores",
                               ylim=(0, float(y_max)+0.2),
                               color='#008ccf' #dark red #8d0034, stjude red #d11947, green #c4d82e, aqua #7ad0e4, dark aqua #17818F
                            )
    
    #remove plot top and left frame border
    #set y label to max @ 1.0
    #add error bars?
    #change color?
    #minor ticks in background

    plt.show()

        
def _main_():
    #find_csv()
    graph_scores()


_main_()