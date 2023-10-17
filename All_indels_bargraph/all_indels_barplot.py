import os
import shutil
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import sys

'''
Changelog: 101623 - fixed rounding error not applying to indel_df.  Line: 101
Changelog: 101723 - removed quick patch from 101623
                  - reordered legend. Line: 147
                  - will plot to 2 decimals but only label whole numbers.  Line: 129

'''



'''
PMH 8/23
Version: 20231017.1
This scripts is to replace the orange bar plots with stacked all indel bar plots.

1) You will first need to run the all_indels version of your results summary.
2) Arrange your samples with WT on top next followed by each guide as they are listed in the summary columns going from left to right.

Name                    Sample              Total   g1  g2  g3
Miller-Plate*Whatever   WT
                        g1
                        g2
                        g3

You can name them whatever you want to make it easier but the program will just label the graphs as WT,g1,g2, etc..
Be sure to get rid of any extra samples that will not need to be charted.  The program will just look for a WT and pool for each guide.

3)Add BP to the end of the all_indels_summary.  You can move it to the barplots_to_be_run folder but the program copies everything there automatically.
  be sure to give each csv a unique name.
4)Follow the prompts.  If you forget to arrange your result summaries before starting the script, there is a forced paused that will prompt you to check your summary csv's
5)The graphs will be generated and labled as whatever the csv name is so again, make sure they have unique names
'''

pd.options.mode.chained_assignment = None  # default='warn'
plot_dir = os.path.join(os.path.dirname(__file__) + "/Barplots_to_be_run/").replace("\\","/")



#This will find any and all files in the NGS folder that have BP at the end of the file.  Even if they are in the original folder
def find_csv():
    
    print(os.getcwd())
    print("Locating files. Please stand by.")
    #will return all the csvs with BP in the name located in the NGS folder and copy to barplot folder
    csv_list = (glob.glob("**/*BP.csv",recursive=True))
  
    
    
    if len(csv_list) >0:
        for csv in csv_list:
            print(csv)
            try:
                shutil.copy(csv,plot_dir+csv.split("\\")[-1])
            except:#file is already in barplots folder
                continue
        print("All files moved to Barplots_to_be_run folder")
        return csv_list 
    else:
        print("\n\nNo CSV's found. Make sure to 'BP' to end of file name.\n\n")
        sys.exit()

#TODO check multiple inputs and returns
#will this work with multiple csvs?
def get_indels():
    
    os.chdir(plot_dir)
    result_csv = glob.glob('*.csv')
    
    print(f"graphing: {result_csv}")

    for csv in result_csv:
        print(csv)
            
        graph_title = csv
        tmp = pd.read_csv(csv)
        print(tmp)
        #remove uncessary columns
        #need to grab guide columns for bar naming purposes.  Sorts through all columns and uses regex to find 'g' + any number + *. Cats Day3. onto front of each guide name
        guide_columns =  [column for column in tmp.columns if re.match(r'g\d*',column)]
        #insert WT into guide_columns[0].
        guide_columns.insert(0,'WT')
        
        indel_df = tmp[['Out-of-frame','In-frame','0bp']]
        
        print(indel_df.head())
        print(guide_columns)
        
        indel_df.insert(0,'Sample',guide_columns)
        
        print(indel_df)
        
        graph_indels(indel_df,graph_title)
        

def graph_indels(df,title):
    
    df.set_index('Sample',inplace=True)
    
    graph_image_name = title.replace(".csv","")
    
    chart_title = f"gRNA Validation via NGS"
    
    ax = df.plot.bar(
        stacked=True,
        color={"0bp":"black","In-frame":"#8D918B","Out-of-frame":"#c10f3a"}
    )
    
    #sets how far below top of element label is placed
    y_offset = -1
    x_offset = 0
    #label each segement
    for bar in ax.patches:
        if bar.get_height() >= 3:
            ax.text(
            #align middle
            bar.get_x() + bar.get_width() / 2 + x_offset,
            #set label above or below top of element
            bar.get_y() + bar.get_height() / 2 + y_offset,
            #fix rounding 
            np.rint(bar.get_height()).astype(int),
            #horizontal alignement
            ha='center',
            color='white',
            weight='bold',
            size=8
    )   
            
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    #rotate xticks
    plt.title(label=chart_title,ha='center')
    plt.ylabel(f"% editing")
    plt.xticks(rotation=0, rotation_mode='anchor',ha='center')
    ax.xaxis.labelpad = 10
    
    
    #set legend to outside of plot area and set order to Oof,If,0bp
    
    leg_order = [2,1,0]
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[idx] for idx in leg_order],[labels[idx] for idx in leg_order],bbox_to_anchor=(1.4,0.5), loc='center right', borderaxespad=0)
    ax.plot()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}\\{graph_image_name}_barplot.png")
    
    print("All graphs have been generated.")
    #plt.show()
    



find_csv()
#input("Ensure all csv's are formatted properly and press enter to continue.")
get_indels()





