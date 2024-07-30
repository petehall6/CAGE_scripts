import os
import subprocess
import pandas as pd

#Change inputs as needed. Be sure to select the sorted (all guides list).  Files must be .xlsx

merge_sort_list_1 = '140_sorted.xlsx'

merge_sort_list_2 = '140_sorted_additl.xlsx'

merged_output = 'lib140_combo_sorted.xlsx'

draft_output = 'lib140_combo_draft.xlsx'

guide_quota = 5

stars = "*"*60
def library_drafter(sorted_file_name, guide_quota, ntc_df,output_file):
      
    #read in xlsx as df
    sorted_df = pd.read_excel(sorted_file_name)
    
    
    picked_list=[]
    ntc_list=[]
    cur_gene=""
    prev_gene=""
    guide_count=0
    
    for row in sorted_df.itertuples():
        cur_gene=row[2]
        
        if prev_gene == cur_gene and cur_gene != "NTC" and guide_count < guide_quota:#adding until quota filled
            picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            prev_gene = cur_gene
            guide_count +=1
            
        elif prev_gene == cur_gene and cur_gene != "NTC" and guide_count >= guide_quota: #quota hit
            None
            
        elif cur_gene == "NTC":#parse out NTC df
            ntc_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            
        else: #first time new gene encountered
            guide_count=1
            picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            prev_gene=cur_gene
            
    
    picked_df = pd.DataFrame(picked_list,columns=['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3'])
    ntc_df = pd.DataFrame(ntc_list,columns=['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3'])
    picked_df.index = picked_df.index+1
    ntc_df.index = ntc_df.index+1
    
    

    picked_df.sort_values(by=['Name'],inplace=True)
    #Creates a series from Name column (gene.g#).  Creates lambda that splits the seires @ '.g', takes the #side and converts that series as int.
    #s stays a series until the very end when its cast as in with astype.  Have to treat the entire thing as a series

    library_df = pd.concat([picked_df,ntc_df],ignore_index=True)
    library_df.index = library_df.index + 1
    
    
    print(stars)
    print(stars)
    print(library_df.head(20))
    print(library_df.tail(20))
    print(stars)
    print(f"Library Draft Complete.  Please see file {output_file}")
    print(stars)
    
    
    library_df.to_excel(output_file,index=False)
    
    return library_df.shape[0]






merge1_df = pd.read_excel(merge_sort_list_1)
merge2_df = pd.read_excel(merge_sort_list_2)


#creates separate ntc df in case they are missing
try:
    ntc_1_df = merge1_df.loc[merge1_df['Name'].str.contains('NTC')]
    ntc_2_df = merge2_df.loc[merge1_df['Name'].str.contains('NTC')]
    ntc_concat_df = pd.concat([ntc_1_df,ntc_2_df])
except:
    ntc_concat_df = ntc_1_df #skip ntc_dfs

merged_concat_df = pd.concat([merge1_df,merge2_df],ignore_index=True)

#sort guides before catting ntc_df
merged_concat_df = merged_concat_df.sort_values(by=['Name','Long_0','Long_1','Long_2','Long_3'],ignore_index=True)

sorted_df = pd.concat([merged_concat_df, ntc_concat_df],ignore_index=True)
sorted_df.index = sorted_df.index + 1

#shouldn't be duplicates but just to be safe
sorted_df.drop_duplicates(inplace=True,keep='last')

sorted_df.to_excel(merged_output, index=False)

print(sorted_df)


#rerun quota for draft
combined_library_size = library_drafter(merged_output,guide_quota,ntc_concat_df,draft_output)

print(f"Initial size number of guides in {merge_sort_list_1}: {merge1_df.shape[0]}")
print(f"Initial size number of guides in {merge_sort_list_2}: {merge2_df.shape[0]}")
print(f"\n\nTotal number of guides in {draft_output}: {combined_library_size}")