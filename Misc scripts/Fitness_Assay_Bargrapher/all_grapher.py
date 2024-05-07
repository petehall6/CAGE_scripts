import os
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt


def label_bar(bars):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2.0, 1.0*height,
                            f'{height}',
                            ha='center', va='bottom')


csv = 'CelFi_ratios_all.csv'

ratio_df = pd.read_csv(csv)

print(ratio_df.head())

fig, ax = plt.subplots()

#cage_name = [c_num for c_num in ratio_df['CAGE']]
gene_name = [g_name for g_name in ratio_df['Gene']]
guide_name = [name for name in ratio_df['Guide']]

print(guide_name)

ratio = [ratio for ratio in ratio_df['Ratio']]

chart_title = "Cell Fitness (CelFi) Assay"

    
bars = ax.bar(gene_name,
                ratio, 
                color='#16A085', 
                edgecolor='#17202A',
                #width=0.5,
                )                  
#label_bar(bars)
#plt.xticks(guide_name)

#plt.ylim(ymax=1.5, ymin=0)
plt.xlabel('Fitness Ratio')
plt.title(chart_title)
#plt.xticks(rotation=10)
ax.spines[['top','right']].set_visible(False)
plt.show()
   