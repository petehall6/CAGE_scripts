import pandas as pd

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


indel_format = wb.add_format({'bg_color':'pink'})

snp_format = wb.add_format({'bg_color': 'red'})

snp_condition = ws.conditional_format('W2:W341',
                      {'type': 'cell',
                       'criteria': '>',
                       'value': '1.03',
                       'format': snp_format})

#Indel columns: G,I,K,M,O,Q,s,U
#=MOD(COLUMN(),3)=0

in_frame = ws.conditional_format('G2:G385',
                                 {'type': 'formula',
                                  'criteria': '=AND(ABS(MOD($G2,3))=0,$G2<>0)',
                                  'format': indel_format})

df.to_excel(writer)

writer.save()

print("Format complete")



