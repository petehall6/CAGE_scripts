import pandas as pd
import openpyxl
hw_results = "STR automation testing-Heartwell output.xlsx"

df_hw_results = pd.read_excel(hw_results, engine='openpyxl')
df_hw_results = df_hw_results.astype(str)
df_hw_results.columns = df_hw_results.columns.str.strip()



#difference returns columns in df NOT listed. So only input columns you want to keep
df_hw_results = df_hw_results.drop(df_hw_results.columns.difference(['Sample Name', 'Marker', 'Allele 1', 'Allele 2', 'Allele 3']), axis=1)
#print(df_hw_results.head(20))
df_hw_results.columns.str.strip(' ')
df_hw_results.replace("Papizan","", regex=True, inplace=True)
df_hw_results.insert(1,'Alleles',df_hw_results[['Allele 1','Allele 2','Allele 3']].agg(','.join,axis=1),True)
df_cleaned = df_hw_results.drop(df_hw_results.columns.difference(['Sample Name','Marker','Alleles']), axis=1)
#print(df_cleaned.head())
df_cleaned['Sample Name'].str.strip()
df_cleaned['Marker'].str.strip()
df_cleaned['Alleles'].str.strip()

print(df_cleaned.head())

print("\nPivoted\n")

df_pivot = df_cleaned.pivot(index='Sample Name', columns='Marker', values='Alleles')

#df_pivot = df_pivot.drop(df_pivot['Alleles'])
print(df_pivot.head())



df_pivot.replace(r' ','',regex=True, inplace=True)
df_pivot.replace(r',$','', inplace=True, regex=True)#regex matches end of string-removes all cases of {allele},
df_pivot.replace(r',$','', inplace=True, regex=True)#removes all cases of {allele},,
df_pivot.replace(' Marker','', inplace=True)
df_pivot = df_pivot.reset_index()

print(df_pivot.head())

df_pivot.columns.str.strip(' ')


df_pivot.to_excel('pivoted_STR_results.xlsx', index=False)