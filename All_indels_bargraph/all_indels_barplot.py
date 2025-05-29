import os
import shutil
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
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
Changelog: 030624 - added main namespace for importing
Changelog: 062824 - guide name columns are now regex'd by search instead of match.  Will resolve any guide names that have the gene/castype in the name
Changelog: 080824 - remove tight layout and added in bbox='tight' to savefig().  Fixed bug when sample names were very long and causing legend to display overtop graph 
Changelog: 030525 - added in a celfi_title variable to change the title of the graph if it is a CelFi.
Changelog: 050125 - added 'ssODN' to column exclusion so ssODN columns will not be needed to be removed when running ssODN barplots.
Changelog: 050925 - fixed bug where hyphen in cell line would cause the x label to be incorrect.                    
Changelog: 052925 - added in fitness graphing function.  Will graph fitness data if 'fitness' is in the file name.

'''



'''
PMH 11/23

Version: 2.0_20250305

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

***
IF running a CelFi include the cell type in ALL CAPS


2)The graphs will be generated and labled as whatever the csv name is so again, make sure they have unique names

3)Completed graphs will be placed in the Barplots_to_be_run folder at the top of the joined folder
'''

#This will find any and all files in the NGS folder that have BP at the end of the file.  Even if they are in the original folder

pd.options.mode.chained_assignment = None
plot_dir = os.path.join(os.path.dirname(__file__) + "/Barplots_to_be_run/").replace("\\","/")


def main():
    find_csv()
    get_indels()

def find_csv():
    
    #print(os.getcwd())
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
    failed_list = []
    fitness = False
    for csv in result_csv:
        try: 
        
            OoF_stand_alone = False    
            graph_title = csv
            tmp = pd.read_csv(csv, encoding="utf-8")
            print(tmp)

            
            #activity's have multiple guides. Sorts through all columns and uses regex to find 'g' + any number
            sample_columns =  [column for column in tmp.columns if re.search(r'\.g\d*|\Ag\d*',column) and "%" not in column and "ssODN" not in column]
           #input(sample_columns)
            
            #insert WT into guide_columns[0].
            sample_columns.insert(0,'WT')
           # input(sample_columns)
            indel_df = tmp[['Out-of-frame','In-frame','0bp']]
            
            #If running for guide activity, the number of guides and samples (not inluding WT) should be 1:1
            #If they are not it will raise an exception
            #If samples and guides are mismatched I'm assuming that they want to plot an edited pool
            #Takes exception, resets the sample names, takes name from the 'Sample' column instead of the guide columns and reinserts
            
            if "activity" in csv.lower():
                    indel_df.insert(0,'Sample',sample_columns)
                    tilt = False
                    celfi_title = False

                    
            #celfi and pools will use sample names
            else:
                sample_columns.clear()
                sample_columns = tmp['Sample'].tolist()
                indel_df.insert(0,'Sample',sample_columns)
                tilt = True
                celfi_title = False
                
            if "celfi" in csv.lower():
                OoF_stand_alone = True
                celfi_title = True

            if "fitness" in csv.lower():
                fitness = True
                
            
            if fitness == True:
                fitness_df = indel_df[['Sample','Out-of-frame']]
                tilt = True
                graph_fitness(fitness_df, graph_title, tilt)
                
            else:
                graph_indels(indel_df, graph_title, tilt, OoF_stand_alone, failed_list, celfi_title)
            os.remove(csv)
            
        except Exception as e:
            print(f"Error with file {csv}: {e}")
            failed_list.append(csv)
            print(f"The following files did not graph: {failed_list}")
            try:
                for csv in failed_list:
                    print(f"File: {csv} did not plot and has been deleted from the Barplots_to_be_run folder")
                    os.remove(csv)
            except:
                pass


def get_cell_line_title(df, gene):
        #finds first row with 'g' in the index and removes the guide and date info
        sample_info = str(df[df.index.str.contains("g")].iloc[0].name).strip()
        cell_line = str(re.sub(r'^d\d*.g\d*|^g\d*.d\d*.|\sg\d*.d\d*.|\sd\d*.g\d*',"", sample_info,flags=re.IGNORECASE)).upper().strip()
        
        #removes cell line from the index and strips whitespace
        df.index = df.index.str.replace(cell_line,"",regex=True,flags=re.IGNORECASE).astype(str)
        df.index = df.index.str.strip().astype(str)
        #add cell_line back into intial row so the x label isnt blank.  Also forces lower case for guide/day 
        index_list = df.index.to_list()
        
        index_list = [index.lower() for index in index_list]
        #sets cell_line as the WT/first sample name
        index_list[0] = cell_line
        
        df.index = index_list
    
        
        chart_title = f"{gene} CelFi assay - {cell_line}"
        
        
        graphing_df = df.copy()
    
        return chart_title, graphing_df

def graph_indels(df, graph_title,tilt, OoF_stand_alone, failed_list,celfi_title):
    
    def _graph_oof():
        df.drop(columns=["0bp","In-frame"], inplace=True)
        #print(f"Dropped columns: {df}")
        graph_image_name = graph_title.replace(".csv","")
        
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

    df.set_index('Sample',inplace=True)
    
    graph_image_name = graph_title.replace(".csv","")
    
    gene = graph_title.split("_")[1]
    print(f"Celfi: {celfi_title}")
    print(f"Gene: {gene}")
    
    if celfi_title == True:
        chart_title, graphing_df = get_cell_line_title(df, gene)

    else:
        chart_title = f"gRNA Validation via NGS"
        graphing_df = df.copy()
        
    ax = graphing_df.plot.bar(
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

    plt.savefig(f"{plot_dir}\\{graph_image_name}_barplot.png",bbox_inches='tight',pad_inches=0.2)
    
    if OoF_stand_alone == True:
        _graph_oof()
    
    
    print("\n\nGraphing complete.")
    
    if len(failed_list) > 0:
        print(f"The following projects failed to graphs. Double check formating. For activity plots, the number of guides must equal the number of samples (not including WT). {failed_list}\n\n")
    else:
        print("Graphing complete.")

def graph_fitness(fitness_df, graph_title, tilt):
    
    def _find_max_y(fitness_score_df):
    
        y_max = fitness_score_df['fitness_score'].max()
            
        if y_max > 1:

            y_buffer = 1.05
            
            y_upper_bound = y_buffer * y_max
            
        else:
            y_upper_bound = 1

        return y_upper_bound
    

    #inital samples will be even number rows, final will be odd number rows
    initial_df = fitness_df.iloc[::2, :]
    final_df = fitness_df.iloc[1::2, :]
    
    #parse out guide names from the sample column
    initial_df['Sample'] = initial_df['Sample'].str.lower().str.extract(r'(g\d*)')
    final_df['Sample'] = final_df['Sample'].str.lower().str.extract(r'(g\d*)')
    
    initial_df = initial_df.rename(columns={'Out-of-frame':'init_oof'})
    final_df = final_df.rename(columns={'Out-of-frame':'final_oof'})
    
    fitness_score_df = pd.merge(initial_df, final_df, on='Sample', how='outer')
    
    fitness_score_df['fitness_score'] = (fitness_score_df['final_oof'] / fitness_score_df['init_oof']).round(2)
    fitness_score_df = fitness_score_df.drop(columns=['init_oof', 'final_oof'])
    
    print(fitness_score_df)
    
    #tries to standardize the bar width.  Bar width is relative to the plot area so the more bars the smaller the width.
    #Fewer bars look weirdly wide.  Set width to 0.3 if 4 or fewer bars looks better, more than five increase
    #width so they don't look too thin
    
    bar_num = len(fitness_score_df['Sample'].to_list())
    y_limit= _find_max_y(fitness_score_df)
    
    if bar_num < 4:
        bar_width = 0.3
    else:
        bar_width = 0.5
    
    fitness_plot = fitness_score_df.plot.bar(
                    x='Sample',
                    y='fitness_score',
                    rot=0,
                    legend=False,
                    ylim=(0, y_limit), 
                    color='#008ccf',
                    align = 'center',
                    label = 'yes',
                    #sets bars in front of grid lines, effect doesn't work unless at least 3
                    zorder = 3, 
                    width = bar_width
                )
    
    #writes labels on top of bars
    for bar in fitness_plot.patches:
        fitness_plot.text(
            #align middle
            bar.get_x() + bar.get_width() / 2,
            #set label slightly above bar
            bar.get_height() + 0.02,
            #label string, and keep 2 decimals
            '{:.2f}'.format(bar.get_height()),
            #horizontal alignement
            ha='center',
            color='black',
            weight='bold',
            size=12
        )

    #turns of the top and right border of the graphing area
    fitness_plot.spines['top'].set_visible(False)
    fitness_plot.spines['right'].set_visible(False)
    
    #draw grid and set behind the bars
    fitness_plot.grid(axis='y', which='both', zorder=3)
    fitness_plot.yaxis.set_minor_locator(AutoMinorLocator(2))

    graph_image_name = graph_title.replace(".csv","")
    project_number = graph_title.split("_")[0]
    gene = graph_title.split("_")[1]
    
    project_title = f"Fitness Scores {project_number} {gene}"

    print(project_title)

    
    #more y axis minor ticks
    plt.xticks(rotation=0, rotation_mode='anchor',ha='center')
    plt.xlabel('Guide',fontsize=15, weight='bold',labelpad=10)
    plt.ylabel('Fitness Score',fontsize=15, weight='bold',labelpad=10)
    plt.title(project_title, fontsize=20, weight='bold', pad=10, ha='center')
    
    plt.savefig(f"{plot_dir}\\{graph_image_name}_barplot.png",bbox_inches='tight',pad_inches=0.2)
    
    print("\n\nGraphing complete.")

if __name__ == "__main__":

    main()