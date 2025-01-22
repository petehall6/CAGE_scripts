import pandas as pd

offinder_df = pd.read_csv('mouse_cas9_offinder.txt',sep='\t') #only has high ota guides


#input(offinder_df)
combined_guides_df = pd.read_csv('mouse_cas9_combo.txt',sep='\t')
#input(combined_guides_df)

tmp_df = pd.concat([offinder_df,combined_guides_df])

true_zero_df = tmp_df.drop_duplicates(subset='sequence', keep=False)

true_zero_df = true_zero_df.drop(columns=['gene','Long_0','Long_1','Long_2','Long_3'])

print(true_zero_df)

true_zero_df.to_csv('true_zero_df_human_cas_9.txt', index=False, sep='\t')