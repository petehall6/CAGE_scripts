import pandas as pd
import os


fpath = "C:\\Users\\phall3\\Documents\\repos\\CAGE_scripts\\NGS_auto_format\\sorted_xl.xlsx"

wb = pd.ExcelFile("raw_xl.xlsx")
df = wb.parse("Sheet1")


#sorting/filtering needs to be done via pandas
#conditional formatting will be easier with openpyxl
writer = pd.ExcelWriter("sorted_xl.xlsx", engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1')
wb = writer.book
ws = writer.sheets['Sheet1']



df = df.sort_values(by="Total_indel",ascending=False)
#explicitly cast total indel column as float for percentages
df['Total_indel'] = df['Total_indel'].astype(float).map('{:.3%}'.format)

#TODO split indel reads into separate columns for easier sorting?

snp_out_format = wb.add_format({'bg_color':'green'})



#out of frame indels in reads 1 and 2 and SNP between 0.98 and 1.03 are highlighted 3
big_conditional = ws.conditional_format('A2:Y341',
                      {'type': 'formula',
                       'criteria': '=AND($W2>0.98,$W2<1.03,ABS(MOD($G2,3))<>0,$G2<>0)',
                       'format': snp_out_format})



df.to_excel(writer)

writer.save()
#for os open
writer.close()
print("Format complete")
os.system(fpath)



