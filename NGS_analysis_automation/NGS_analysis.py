import pandas as pd
import os

#Create PCR Analysis folder with outside PCR request

ngs_df = pd.read_excel('NGS_SRM.xls')

columns_drop = ['Date Ordered','Gene Name/Gene ID','Number of Tube Samples/Plates','Sample Format','Specify Sample Type','SRM Sample #','Consolidation Plate?','Number of Consolidation Plates','User Comments','Lab Comments']

ngs_df.drop(columns_drop, axis=1,inplace=True) 

outside_df = ngs_df[ngs_df['Principal Investigator'] != 'Miller, Shondra']
outside_df = outside_df[ngs_df['Sample Type'] !='Tail Snip/Toe Snip']


print(outside_df)

os.chdir('PCR_Analysis')

for idx, row in outside_df.iterrows():
    
    dir_name = str(outside_df['CAGE Project #'][idx]) + '_' + str(outside_df['SRM Order #'][idx])
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

