import pandas as pd

#File said corrupted.  When I tried to open with pd.read_excel(file, engine=openpyxl)
#I got an error saying something like cant open Openoffice ODS
#installed pyODF library and ran this
#can get conda to install openoffice library

#
file1 = "Suppl_TableS1.xlsx"
file2 = "Suppl_TableS2.xlsx"
file3 = "Suppl_TableS3.xlsx"
file4 = "Suppl_TableS4.xlsx"
file5 = "Suppl_TableS5.xlsx"



df1 = pd.read_excel(file1, engine="odf")
df2 = pd.read_excel(file2, engine="odf")
df3 = pd.read_excel(file3, engine="odf")
df4 = pd.read_excel(file4, engine="odf")
df5 = pd.read_excel(file5, engine="odf")

df1.to_excel("output_table1.xlsx")
df2.to_excel("output_table2.xlsx")
df3.to_excel("output_table3.xlsx")
df4.to_excel("output_table4.xlsx")
df5.to_excel("output_table5.xlsx")

print("Completed")