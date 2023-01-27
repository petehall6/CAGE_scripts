import pandas as pd

wb = pd.ExcelFile("raw_xl.xlsx")
df = wb.parse("Sheet1")

def sort_indel(df):
    df = df.sort_values(by="Total_indel",ascending=False)
    #explicitly cast total indel column as float for percentages
    df['Total_indel'] = df['Total_indel'].astype(float).map('{:.3%}'.format)
    writer = pd.ExcelWriter("sorted_xl.xlsx")
    df.to_excel(writer)
    writer.save()






                    
df.style.apply(lambda x: ["background-color: red" 
                          if (i == 23 and (v > 1.03)) 
                          else "background-color: blue" for i, v in enumerate(x)], axis = 1)
                    
                    
                    

writer = pd.ExcelWriter("sorted_xl.xlsx", engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1')
wb = writer.book
ws = writer.sheets['Sheet1']
snp_format = wb.add_format({'bg_color': 'red'})

ws.conditional_format('W2:W341',
                      {'type': 'cell',
                       'criteria': '>',
                       'value': '1.03',
                       'format': snp_format})

writer.save()



