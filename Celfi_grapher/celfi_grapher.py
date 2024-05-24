import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil
import glob
import re

inputCSV = "input.csv"
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
script_dir = os.getcwd()
plt.rcParams["savefig.directory"] = (os.path.join(os.environ["USERPROFILE"], "Desktop"))

csv_dir = os.path.join(os.getcwd(),'csv files')
print(csv_dir)

def find_csv():
    #clear csv_dir
    def clear_csv_dir():
        os.chdir(csv_dir)
        
        for file in os.scandir(csv_dir):
            os.remove(file)
        os.chdir(script_dir)
    
    clear_csv_dir()
    #get input from csv
    input_df = pd.read_csv(inputCSV)

    input_df = input_df.astype(str)
    
   #input(input_df.head())
    
    #puts columns into list to containerize the data eaiser
    input_projects = input_df.values.tolist()
    print("Scanning for csv")
    
    for target_proj in input_projects:
        
        cage_proj, ngs_date = target_proj
        
        #appends 0 in from of months Jan->Sept
        if len(ngs_date) < 6:
            ngs_date = '0' + str(ngs_date)
        
        
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined').replace("\\\\","\\")
        os.chdir(joined_dir)
        
        #search joined folder for all_indel
        for file in glob.glob(f"{cage_proj}**\*all_indels.csv",recursive=True):
            if cage_proj in file:
                shutil.copy(file, csv_dir)
    
    print("CSV files copied over.")

def get_scores():
    def _pick_guides(comp):            
            #*debug
            '''
            if comp ==1 :
                initial_time_pt_choice = 1
                final_time_pt_choice = 2
            else:
                initial_time_pt_choice = 3
                final_time_pt_choice = 5
            '''
            
            #print(f"Comparions #{comp}: ")
            initial_time_pt_choice = int(input(f"initial time point for comparison {comp}: "))
            final_time_pt_choice = int(input(f"final time point for comparison {comp}: "))
            
            #gives list [sample, g, guide_num]
            init_time_name = oof_df.loc[initial_time_pt_choice][0]
            init_time_guide = re.findall('g\d*',init_time_name)[0]
            init_time_oof = (oof_df.loc[initial_time_pt_choice][1])
            
            final_time_name = oof_df.loc[final_time_pt_choice][0]
            final_time_guide = re.findall('g\d*',final_time_name)[0]
            final_time_oof = (oof_df.loc[final_time_pt_choice][1])

            #checks to see if guides match
            
            if init_time_guide != final_time_guide:
                input("Guide names don't match.  Press Enter to repick guides.")
                while init_time_guide != final_time_guide:
                    initial_time_pt_choice = int(input(f"initial time point for comparison {comp}: "))
                    final_time_pt_choice = int(input(f"final time point for comparison {comp}: "))
                
                    #gives list [sample, g, guide_num]
                    init_time_name = oof_df.loc[initial_time_pt_choice][0]
                    init_time_guide = re.match('g\d*',init_time_name).group(0)
                    init_time_oof = (oof_df.loc[initial_time_pt_choice][1])
                    
                    final_time_name = oof_df.loc[final_time_pt_choice][0]
                    final_time_guide = re.match('g\d*',final_time_name).group(0)
                    final_time_oof = (oof_df.loc[final_time_pt_choice][1])

            guide_name = init_time_guide
            fa_score = round(final_time_oof/init_time_oof,2)
            
            
            return guide_name, fa_score
        
    fa_score_list=[]
    
    print(f"this is the fa-score list: {fa_score_list}")
    
    os.chdir(csv_dir)
    
    for csv in os.scandir(csv_dir):
        proj_name = csv.name.split("_")[0]
        gene_name = csv.name.split("_")[1]
        print(proj_name)

        oof_columns = ['Sample','Out-of-frame']
        oof_df = pd.read_csv(csv,usecols=oof_columns)
        oof_df.index = oof_df.index+1

        print(oof_df)
        
        #*debug
        #comp_num=3
        
        comp_num = int(input("Number of Comparisons: "))
        
        for comp in range(1, comp_num+1):
            
            guide_name, fa_score = _pick_guides(comp)
            
            print(f"\n\nFitness assay score for comparison {comp}: {fa_score}\n\n")
            
            fa_score_list.append([guide_name,fa_score])
            
        print(f"Fitness scores: {fa_score_list}")
        
        return fa_score_list, proj_name, gene_name

def graph_scores():

    fa_scores, proj, gene = get_scores()
    
    #convert to dataframe to access df.plot.bar(), gives more customization
    fa_score_df=pd.DataFrame(fa_scores, columns=["guide","fa_score"])
    
    bar_num = len(fa_score_df['guide'].to_list())
    
    #tries to standardize the bar width.  Bar width is relative to the plot area so the more bars the smaller the width.
    #Fewer bars look weirdly wide.  Set width to 0.3 if 4 or fewer bars looks better, more than five increase
    #width so they don't look too thin
    if bar_num < 4:
        bar_width = 0.3
    else:
        bar_width = 0.5

    print(fa_score_df)
    
    fa_plot = fa_score_df.plot.bar(
                     x="guide",
                     xlabel = "sgRNA",
                     y='fa_score',
                     ylabel = 'Fitness Score',
                     rot=0,
                     legend=False,
                     title=f"{proj}.{gene}",
                     ylim=(0, 1.1), 
                     color='#008ccf',
                     align = 'center',
                     label = 'yes',
                     #sets bars in front of grid lines, effect doesn't work unless at least 3
                     zorder = 3, 
                     width = bar_width
                )
    
    #writes labels on top of bars
    for bar in fa_plot.patches:
        fa_plot.text(
            #align middle
            bar.get_x() + bar.get_width() / 2,
            #set label slightly above bar
            bar.get_height() + 0.02,
            #label string, and keep 2 decimals
            "{:.2f}".format(bar.get_height()),
            #horizontal alignement
            ha='center',
            color='black',
            weight='bold',
            size=8
        )

    #turns of the top and right border of the graphing area
    fa_plot.spines['top'].set_visible(False)
    fa_plot.spines['right'].set_visible(False)
    
    #draw grid and set behind the bars
    fa_plot.grid(axis='y', zorder=0)

    plt.show()


def main():
    find_csv()
    graph_scores()

main()