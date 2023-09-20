import os
from pathlib import Path
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


#set saving directory
graph_dir = "Z:ResearchHome\\Groups\\millergrp\\home\\common\\NGS\\082223\\joined\\_cell_fitness_folder\\good_graphs\\even_better_graphs"

mpl.rcParams["savefig.directory"] = Path(graph_dir)


def label_bar(bars):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2.0, 1.0*height,
                            f'{height}',
                            ha='center', va='bottom')



csv = '_all_ratios_x6.csv'
chart_title = "Cell Fitness (CelFi) Assay - HGGx6"
input(print(os.getcwd()))
os.chdir('csv_files/ratios/')



ratio_df = pd.read_csv(csv)

print(ratio_df.head())

fig, ax = plt.subplots()

#cage_name = [c_num for c_num in ratio_df['CAGE']]
gene_name = [g_name for g_name in ratio_df['Gene']]
guide_name = [name for name in ratio_df['Guide']]

print(guide_name)

ratio = [ratio for ratio in ratio_df['Ratio']]



what_we_are_plotting = gene_name
    
bars = ax.bar(what_we_are_plotting,
                ratio, 
                color='#16A085', 
                edgecolor='#17202A',
                #width=0.5,
                )                  
#label_bar(bars)
plt.xticks(what_we_are_plotting,rotation=90)
plt.tight_layout()
plt.ylim(ymax=1.6, ymin=0)
plt.xlabel('Fitness Ratio')
plt.title(chart_title)
#plt.xticks(rotation=10)
ax.spines[['top','right']].set_visible(False)
plt.show()
   