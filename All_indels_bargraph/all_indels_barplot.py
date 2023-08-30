import os
import matplotlib.pyplot as plt
import pandas as pd
import re
import sys
import glob


#Replacing orange bar plots






#right now just finds all csv's with BP in name in current directory.  Will need to adapt to run on entire NGS folder and BP to run folder
def find_csv(path=os.getcwd()):

    os.chdir('test_Barplots_to_be_run/')
    csv_list=glob.glob('*.csv')
    
    #will return all the csvs with BP in the name located in the NGS folder
    csv_list = glob.glob("*/*testBP.csv")

    print(csv_list)
    input()
    if len(csv_list) >0:
        return csv_list 
    else:
        print("\n\nNo CSV's found. Make sure to 'BP' to end of file name.\n\n")
        sys.exit()

#TODO check multiple inputs and returns
#will this work with multiple csvs?
def get_indels(result_csv):
    for csv in result_csv:
        
        graph_title = csv
        tmp = pd.read_csv(csv)
        
        #remove uncessary columns
        #need to grab guide columns for bar naming purposes.  Sorts through all columns and uses regex to find 'g' + any number + *. Cats Day3. onto front of each guide name
        guide_columns =  ['Day3.'+column for column in tmp.columns if re.match(r'g\d*',column)]
        #insert WT into guide_columns[0].
        guide_columns.insert(0,'WT')
        
        indel_df = tmp[['Out-of-frame','In-frame','0bp']]
        indel_df.insert(0,'Sample',guide_columns)
        #convert indel columns to int to round to whole number
        indel_df = indel_df.astype({'Out-of-frame': int, 'In-frame': int, '0bp': int})
        
        print(indel_df)
        
        graph_indels(indel_df,graph_title)
        

def graph_indels(df,title):
    
    df.set_index('Sample',inplace=True)
    
    indel_title = title.split("_")[0]
    
    chart_title = f"gRNA Validation via NGS"
    
    ax = df.plot.bar(
        stacked=True,
        color={"Out-of-frame":"#c10f3a","In-frame":"#8D918B","0bp":"black"},
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
            bar.get_height().astype(int),
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
    plt.xticks(rotation=45, rotation_mode='anchor',ha='right')
    ax.xaxis.labelpad = 10
    #set legend to outside of plot area
    plt.legend(bbox_to_anchor=(1.4,0.5), loc='center right', borderaxespad=0)
    ax.plot()
    plt.tight_layout()
    plt.savefig(f"graphs\\{indel_title}_barplot.png")
    plt.show()
    



csv = find_csv()

get_indels(csv)





