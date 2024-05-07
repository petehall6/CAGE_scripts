import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

"""
    This script regraphs previously graphed data but this time it will lock the y-axis @ 1.2
    Drop all the .CSV's into this folder and run the script
    
"""


def label_bar(bars):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2.0, 1.0*height,
                            f'{height}',
                            ha='center', va='bottom')



dir = os.getcwd()

print(dir)
csv_list=[]

for file in os.scandir(dir):
    if file.name.endswith('.csv'):
        csv_list.append(file)


#print(csv_list)        
for csv in csv_list:
    ratio_df = pd.read_csv(csv)
    
    gene_name = csv.name.split("_")[0]
    guide_name = [name for name in ratio_df['Guide']]
    ratio = [ratio for ratio in ratio_df['Ratio']]

    chart_title = gene_name +" Cell Fitness (CelFi) Assay"
    num_of_guides = len(guide_name)

    if num_of_guides == 1:
        fig, ax = plt.subplots()
        #makes the individual guide bar plot a separate entity for easier editing
        
        x_pos = [0,1,2]
        
        bars = ax.bar(x_pos,
                        [0,ratio[0],0], 
                        color='#16A085', 
                        edgecolor='#17202A',
                        width=0.5,
                        )                  
        label_bar(bars)
        plt.xticks(range(1,2),guide_name)

        plt.ylim(ymax=1.2, ymin=0)
        plt.ylabel('Fitness Ratio')
        plt.title(chart_title, y=1.05)

        ax.spines[['top','right']].set_visible(False)
        plt.show()
    if num_of_guides == 2:
        fig, ax = plt.subplots()
        #makes the individual guide bar plot a separate entity for easier editing
        
        x_pos=[0,1,2,3]
        
        bars = ax.bar(x_pos,
                        [0,ratio[0],ratio[1],0],
                        color='#16A085', 
                        edgecolor='#17202A',
                        width=0.5,
                        )                  
        label_bar(bars)
        plt.xticks(range(1,3), guide_name)

        plt.ylim(ymax=1.2, ymin=0)
        plt.ylabel('Fitness Ratio')
        plt.title(chart_title, y=1.05)
        ax.spines[['top','right']].set_visible(False)

        plt.show()