import os
import shutil
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import re
import sys

'''
Changelog: 101623 - fixed rounding error not applying to indel_df. 
Changelog: 101723 - removed quick patch from 101623
                  - reordered legend. Line: 147
                  - will plot to 2 decimals but only label whole numbers. 
Changelog: 102023 - Default setting is to run guide activity plot and names bars after guides
                  - Added a try/except loop that will handle insertion length error and change
                    the bar names to the 'Sample' column from the tmp df.
Changelog: 111023 - Added tags to file names to better distribute work flow
                  - Running celfi's will also generate an out of frame graph
Changelog: 012224 - Added try/except loop to catch graphs that fail so downstream graphs are not effected. Line: 97
                    


'''



'''
PMH 11/23

Version: 2.0_20231110

This scripts is to replace the orange bar plots with stacked all indel bar plots.

1) You will first need to run the all_indels version of your results summary.

****If you are running it for guide activity/orange bar plot

    add activity BP to the end of the file name
    This program uses this to name the samples after the guides and format the graphs correctly

    CAGE1234 all_indels activity BP.csv

****If you are running an editied cell pool or celfi

    Arrange your samples with WT on top next followed by each guide try to keep the names short so they dont aren't crazy long in the image
    Name                    Sample              Total   g1
    Miller-Plate*Whatever   WT
                            KO pool BT16
                            dKO Pool BT16/Sirt

    add pool  or celfi BP to the end of the file name
    This program uses this to name the samples after the 'Sample' column in the csv and format the graphs correctly

    CAGE1234 all_indels pool BP.csv
    CAGE1234 all_indels celfi BP.csv

2)The graphs will be generated and labled as whatever the csv name is so again, make sure they have unique names

3)Completed graphs will be placed in the Barplots_to_be_run folder at the top of the joined folder
'''

pd.options.mode.chained_assignment = None
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

def get_indels():
    
    os.chdir(plot_dir)
    result_csv = glob.glob('*.csv')
    
    print(f"graphing: {result_csv}")

    for csv in result_csv:
        try: 
        
            OoF_stand_alone = False    
            graph_title = csv
            tmp = pd.read_csv(csv)
            print(tmp)
        
            #activity's have multiple guides. Sorts through all columns and uses regex to find 'g' + any number
            sample_columns =  [column for column in tmp.columns if re.match(r'g\d*',column)]
            #insert WT into guide_columns[0].
            sample_columns.insert(0,'WT')
            indel_df = tmp[['Out-of-frame','In-frame','0bp']]
            
            #If running for guide activity, the number of guides and samples (not inluding WT) should be 1:1
            #If they are not it will raise an exception
            #If samples and guides are mismatched I'm assuming that they want to plot an edited pool
            #Takes exception, resets the sample names, takes name from the 'Sample' column instead of the guide columns and reinserts
            
            if "activity" in csv:
                    indel_df.insert(0,'Sample',sample_columns)
                    tilt = False
                    
            #celfi and pools will use sample names
            else:
                sample_columns.clear()
                sample_columns = tmp['Sample'].tolist()
                indel_df.insert(0,'Sample',sample_columns)
                tilt = True
                
            if "celfi" in csv:
                OoF_stand_alone = True
                
            
            print(indel_df)
            
            graph_indels(indel_df,graph_title, tilt, OoF_stand_alone)
            os.remove(csv)
        except:
            failed_list = []
            failed_list.append(csv)
            print(f"The following files did not graph: {failed_list}")
            for csv in failed_list:
                os.remove(csv)
                         
def graph_indels(df,title,tilt, OoF_stand_alone):
    
    def graph_oof():
        df.drop(columns=["0bp","In-frame"], inplace=True)
        print(f"Dropped columns: {df}")
        graph_image_name = title.replace(".csv","")
        
        chart_title = f"gRNA Validation via NGS"
        
        ax = df.plot.bar(
            stacked=True,
            color={"Out-of-frame":"#c10f3a"}
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
        if tilt == True:
            plt.xticks(rotation=25, rotation_mode='anchor',ha='right')
        else:
            plt.xticks(rotation=0, rotation_mode='anchor',ha='center')
        ax.xaxis.labelpad = 10
        
        
        #set legend to outside of plot area and set order to Oof,If,0bp
        
        plt.legend(bbox_to_anchor=(1.4,0.5), loc='center right', borderaxespad=0)
        ax.plot()
        ax.set_ylim([0,100])
        plt.margins()
        plt.tight_layout()
        plt.savefig(f"{plot_dir}\\{graph_image_name}_OoF_barplot.png")
        
        return


    #*first pass
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
    if tilt == True:
        plt.xticks(rotation=25, rotation_mode='anchor',ha='right')
    else:
        plt.xticks(rotation=0, rotation_mode='anchor',ha='center')
    ax.xaxis.labelpad = 10
    
    
    #set legend to outside of plot area and set order to Oof,If,0bp
    
    leg_order = [2,1,0]
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend([handles[idx] for idx in leg_order],[labels[idx] for idx in leg_order],bbox_to_anchor=(1.4,0.5), loc='center right', borderaxespad=0)
    ax.plot()
    plt.tight_layout()
    plt.savefig(f"{plot_dir}\\{graph_image_name}_barplot.png")
    
    if OoF_stand_alone == True:
        graph_oof()
    
    
    print("All graphs have been generated.")
    #plt.show()


find_csv()
#input("Ensure all csv's are formatted properly and press enter to continue.")
get_indels()





