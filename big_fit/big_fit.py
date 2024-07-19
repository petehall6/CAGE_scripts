import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import os
import shutil
import glob


import xlwings as xw
import textwrap

inputCSV = "input.csv"
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
script_dir = os.getcwd()
plt.rcParams["savefig.directory"] = (os.path.join(os.environ["USERPROFILE"], "Desktop"))
holding_dir = os.path.join(os.getcwd(),'holding')
compiled_excel = "compiled_data.xlsx"

project_title=""

def find_csv():
    #clear holding_dir
    def clear_holding_dir():
        os.chdir(holding_dir)
        
        for file in os.scandir(holding_dir):
            os.remove(file)
        os.chdir(script_dir)
    
    clear_holding_dir()
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
        for file in glob.glob(f"{cage_proj}**all_inde*\*all_indel*.csv",recursive=True):
            if cage_proj in file:
                shutil.copy(file, holding_dir)
    
    print("CSV files copied over.")

def compile_to_excel():
    print("Compiling data into Excel.")
    
    def _format_excel(compiled_df):
                    
            os.chdir(script_dir)
            current_time = datetime.strftime(datetime.today(),"%m%d%Y_%H%M")
            
            
            #compiled_df.reset_index(drop=True, inplace=True)
            
            #*Transfer compiled df to excel
            
            #load workbook
            template_excel = 'template.xlsx'
            app = xw.App(visible=False)
            wb = xw.Book(template_excel)
            ws1 = wb.sheets['Sheet1']
            ws2 = wb.sheets['dont_touch']
            #update workbook
            ws1.range('B2').options(index=False,header=False).value = compiled_df
            
            #*Transfer comparison groups to excel
            #get unique CAGE#s
            
            time_point_list=[]
            
            cage_nums = sorted(set(compiled_df['CAGE#'].tolist()))
            
            for num in cage_nums:
                time_point_list.append(f'{num}_init_tp')
                time_point_list.append(f'{num}_final_tp')
                
            time_point_df = pd.DataFrame(time_point_list)
            
            ws2.range('A2').options(index=False,header=False).value = time_point_df
            
            #close and save as updated_form
            wb.save(compiled_excel)
            wb.close()
            app.quit()
            
    
    tmp_df_list=[]
    
    #instainate big df will add Sample Name column after comipling all csvs
    compiled_columns = ['CAGE#','Well#','0bp','In-frame', 'Out-of-Frame','Comparison_Group']
    compiled_df = pd.DataFrame(columns=[compiled_columns])
    
    os.chdir(holding_dir)
    
    #will return dir entry so can access name via attribute
    csv_list = os.scandir()
    
    #parse data from each csv
    for csv in csv_list:
        
        tmp_columns = ['Name','0bp','In-frame','Out-of-frame']
        
        #copy only the columns needed and append to list to concat into compiled_df later
        tmp_df = pd.read_csv(csv)
        
        print(csv.name)
        
        
        tmp_df = tmp_df[tmp_columns]
        tmp_df.rename(columns={"Name":"Well#"},inplace=True)
        
        cage_num = csv.name.split("_")[0]
        gene = csv.name.split("_")[1]
        tmp_df.insert(0,"CAGE#",cage_num,True)
        tmp_df.insert(2,"Gene",gene,True)
        
        tmp_df_list.append(tmp_df)
        
    #compile all csv df's and reset index, save as excel
    compiled_df = pd.concat(tmp_df_list,ignore_index=True)
    compiled_df.index = compiled_df.index +1
    #print(compiled_df.head(20))
    

    #compiled_df.to_excel(f"compiled_data_{current_time}.xlsx")
    
    _format_excel(compiled_df)

def get_scores():
    
    print("\nChoose comparison groups in the compiled_data.xlsx.\n")
    #TODO uncomment
    #input("\n\nPress Enter to continue\n\n")
    
    #Read in compiled excel, drop rows without a comparison selected and reset index for readability
    columns = ['CAGE#', 'Gene', 'Out-of-frame', 'Comparison_Group']
    score_df = pd.read_excel(compiled_excel, usecols=columns)
    score_df = score_df.dropna(axis=0, how='any')
    
    score_df.reset_index(drop=True,inplace=True)
    score_df.index = score_df.index + 1
    
    #generate fitness scores
    #break out each initial and final time point into seperate df's then merge to have 1 flat combined df.
    init_df = score_df[["CAGE#","Gene","Out-of-frame"]][score_df["Comparison_Group"].str.contains("init")]
    init_df.rename(columns={"Out-of-frame":"init_oof"},inplace=True) #need to have unique column for merge

    final_df = score_df[["CAGE#","Gene","Out-of-frame"]][score_df["Comparison_Group"].str.contains("final")]
    final_df.rename(columns={"Out-of-frame":"final_oof"},inplace=True)    
    combo_df = init_df.merge(final_df,on=["CAGE#","Gene"])
    
    combo_df.insert(4,"fitness_score",(combo_df["final_oof"] / combo_df["init_oof"]).round(2))
    combo_df.drop(columns=["init_oof","final_oof"],inplace=True)
    
    combo_df.index = combo_df.index + 1
    
    combo_df.sort_values(by=['fitness_score'],inplace=True)
    
    combo_df.to_excel("fitness_scores.xlsx")
    
    graph_scores(combo_df)

def graph_scores(combo_df):
    
    def _rotate_labels(plot):
        labels = []
        for label in plot.get_xticklabels():
            text = label.get_text()
            labels.append(text)
            fa_plot.xaxis.set_tick_params(which='major', pad=2)
        plot.xaxis.set_tick_params(which='both', pad=0)
        plot.set_xticklabels(labels, rotation=45, size=10.0, ha='right', rotation_mode='anchor')
    
    def _wrap_labels(plot, width, break_long_words=True):
        
        labels = []
        for label in plot.get_xticklabels():
            text = label.get_text()
            labels.append(textwrap.fill(text, width=width,
                        break_long_words=break_long_words))
        plot.set_xticklabels(labels, rotation=45, size=10.0, ha='right')
    
    def _find_max_y(fa_score_df):
        
        y_max = fa_score_df['fitness_score'].max()
        
        y_buffer = 1.05
        
        y_upper_bound = y_buffer * y_max
        
        return y_upper_bound
    
    fa_score_df = combo_df[['Gene','fitness_score']]
    
    y_limit= _find_max_y(fa_score_df)
    
    bar_num = len(fa_score_df['Gene'].to_list())
    
    #tries to standardize the bar width.  Bar width is relative to the plot area so the more bars the smaller the width.
    #Fewer bars look weirdly wide.  Set width to 0.3 if 4 or fewer bars looks better, more than five increase
    #width so they don't look too thin
    if bar_num < 4:
        bar_width = 0.3
    else:
        bar_width = 0.5
    fa_plot = fa_score_df.plot.bar(
                     x="Gene",
                     y='fitness_score',
                     rot=0,
                     legend=False,
                     #ylim=(0, y_limit), 
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
            size=12
        )

    #turns of the top and right border of the graphing area
    fa_plot.spines['top'].set_visible(False)
    fa_plot.spines['right'].set_visible(False)
    
    #draw grid and set behind the bars
    fa_plot.grid(axis='y', which='both', zorder=3)
    fa_plot.yaxis.set_minor_locator(AutoMinorLocator(2))
    #fa_plot.minorticks_on()
    
    
    
    #_wrap_labels(fa_plot,width=5)
    _rotate_labels(fa_plot)

    #more y axis minor ticks
    
    plt.xlabel("Gene",fontsize=15, weight='bold',labelpad=10)
    plt.ylabel("Fitness Score",fontsize=15, weight='bold',labelpad=10)
    plt.title(f"Fitness Scores {project_title}", fontsize=20, weight='bold', pad=10)
    plt.autoscale()
    plt.show()
    
    

def main():
    #find_csv()
    #compile_to_excel()
    get_scores()

main()

print("Graphing Complete.")