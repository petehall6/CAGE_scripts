from tkinter import filedialog
import os
import pandas as pd
import numpy as np
import itertools
DESKTOP = os.path.join(os.path.join(os.environ['USERPROFILE'], 'Desktop'))
_PETE = 'Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_boostrap_emailer'

#functions
def open_file():
    
    fileytpes = (
        ('Excel files', '*.xls'),
        ('All files', '*.*')
    )
    
    srm_out = filedialog.askopenfilename(
        title = "Select file",
        initialdir= _PETE,
        filetypes=fileytpes  
    )
    
    return srm_out
    
    #try replace using os.[USERPROFILE]+DESKTOP

def df_from_template(template):
       info = pd.read_excel(template)
       
       tmp = pd.DataFrame(info)

       srm_df = tmp[['SRM Order #',
                     'PI',
                     'Requested By',
                     'Project Number',
                     'Project Scope',
                     'Cell Line of Choice',
                     'Project Objective',
                     'Target Gene Name'
                     ]]
       
       
       results = srm_df.values.tolist()
       
       
       print(results)
       input("hold for input")
       
       
       return results