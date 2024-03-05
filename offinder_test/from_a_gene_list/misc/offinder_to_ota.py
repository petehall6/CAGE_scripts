import pandas as pd
import os

#TODO why doesnt this run on the cluster?

#* Change these as needed.  Don't touch anything else.
# run on local machine

os.chdir("input")
print(os.getcwd())
input_file = 'input_test.txt'
#output_file = 'ota_test_1019.csv'
long_mismatch_num = 3


#summary frame serves as the template for the completed mismatch dataframe
#Count_alignment for Cas 9Long
tmp_df =pd.read_csv(input_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None)
tmp_df.drop(['chrom', "strand","name", 'result'],axis=1,inplace=True)
count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()

count_df.columns = count_df.columns.droplevel(1)
count_df.reset_index(0)

count_df.fillna(0, inplace=True)

count_df.columns = ['Long_0','Long_1','Long_2','Long_3']

print(count_df)


count_df['Long_1'] = count_df.iloc[:,0:1].sum(axis=1)
count_df['Long_2'] = count_df.iloc[:,0:2].sum(axis=1)
count_df['Long_3'] = count_df.iloc[:,0:3].sum(axis=1)

print("*"*20)
print(count_df)
print("*"*50)
os.chdir("../")
count_df.to_csv("output_check.csv")