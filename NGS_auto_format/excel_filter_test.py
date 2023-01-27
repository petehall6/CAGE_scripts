import pandas as pd

wb = pd.ExcelFile("exl_filter.xlsx")
df = wb.parse("Sheet1")


df = df.sort_values(by="Total_indel",ascending=False)
#explicitly cast total indel column as float for percentages
df['Total_indel'] = df['Total_indel'].astype(float).map('{:.3%}'.format)


df.style.apply(lambda x: ["background: red" if v > 1.03 else "" for v in df['SNP_test']], axis = 1)
print()
writer = pd.ExcelWriter("sorted_xl.xlsx")
df.to_excel(writer)
writer.save()