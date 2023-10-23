import pandas as pd
#from icecream import ic
def df_from_template(template):
       info = pd.read_excel(template)
       
       tmp = pd.DataFrame(info)

       srm_df = tmp[['PI',
                     'SRM Order #',
                     'Requested By',
                     'Project Number',
                     'Project Scope',
                     'Cell Line of Choice',
                     'Project Objective',
                     'Target Gene Name',
                     'Species',
                     'Comments',
                     'Is this a human pluripotent stem cell (hESC or hiPSC) project?'
        ]]
       
       #results = srm_df.to_numpy().tolist()
       
       #print(srm_df.columns)

       return srm_df

template = "_CAGEServices_Excel Export_shondra.xls"

df = df_from_template(template)



pi_list = list(set(df["PI"].values.tolist()))

print(pi_list)



