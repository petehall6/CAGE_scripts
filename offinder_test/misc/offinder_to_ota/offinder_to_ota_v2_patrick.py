import pandas as pd
import os

#This version uses grouping instead of itterows. Works much faster

os.chdir("input")
print(os.getcwd())
input_file = 'input_test.txt'
long_mismatch_num = 3
#summary frame serves as the template for the completed mismatch dataframe
#Count_alignment for Cas 9Long


df =pd.read_csv(input_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None)
#df=df.drop_duplicates(subset='gRNA', keep='first')
df.sort_values(by=['gRNA', 'mismatch'])

counted_df = df[['gRNA','pos','mismatch',]].groupby(by=['gRNA', 'mismatch'], axis=0).count().unstack()
counted_df=counted_df.add_prefix('Long_').reset_index()

print(counted_df)


counted_df.to_csv('OTA_scores.csv', index=False) 

counted_df=pd.read_csv('OTA_scores.csv', skiprows=2, names=['gRNA', 'Long_0', 'Long_1','Long_2','Long_3'])
counted_df=counted_df.fillna(0)
print(counted_df.head(8))
counted_df=pd.merge(counted_df, df[['gRNA','name']].drop_duplicates(), on='gRNA', how='left')
counted_df['Long_3']=counted_df['Long_3']+counted_df['Long_2']+counted_df['Long_1']+counted_df['Long_0']
counted_df['Long_2']=counted_df['Long_2']+counted_df['Long_1']+counted_df['Long_0']
counted_df['Long_1']=counted_df['Long_1']+counted_df['Long_0']
print('new counted')
print(counted_df.head(8))
counted_df['Snapgene_label']=counted_df['name']+'_'+counted_df['Long_0'].astype('Int64').astype('str')+'_'+counted_df['Long_1'].astype('Int64').astype('str') \
        +'_'+counted_df['Long_2'].astype('Int64').astype('str')+'_'+counted_df['Long_3'].astype('Int64').astype('str')
counted_df.to_csv('OTA_scores_final.csv', index=False)
