import os
import matplotlib as plt
import numpy as np
import pandas as pd
import glob
import re

def find_csv(path=os.getcwd()):

    csv_list = [csv for csv in os.listdir(path) if csv.endswith('.csv')]
    print(csv_list)
    return csv_list

def remove_parenthetical(dirty_df):
    clean_df = pd.DataFrame()
    
    for column in dirty_df:
        clean_df[column] = dirty_df[column].replace(r"\(.*\)",'',regex=True)
    
    return clean_df

def get_indels(result_csv):
    
    for csv in result_csv:
        
        indel_df = pd.read_csv(csv)
        #drop unneeded columsn
        indel_df.drop(['Name','Total','Total_indel','SNP_test','raw_wt_counter','0bp','In-frame','Out-of-frame'], inplace=True, axis=1)
        #drop last column
        indel_df.drop(indel_df.columns[-1], inplace=True, axis=1)
        indel_df = indel_df.fillna('-')
        print(indel_df.head())
 
        scrubbed_df = remove_parenthetical(indel_df)
        
        print(scrubbed_df)
        
        #find all the indel columns
        indel_columns = [column for column in scrubbed_df.columns.values.tolist() if 'Indel' in column]
        
        #loop through list of columns to find oof, if, and 0bp indels
        
        
        
        
        
        
        print(indel_columns)
        
        
        
        
        
        
        
        
        
        
        
    return


result_csv = find_csv()

get_indels(result_csv)
