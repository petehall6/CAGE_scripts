from tkinter import filedialog
from tkinter import ttk
import os
import pandas as pd
import numpy as np
import itertools
import re

DESKTOP = os.path.join(os.path.join(os.environ['USERPROFILE'], 'Downloads'))
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
        #initialdir= DESKTOP,
        filetypes=fileytpes  
    )
    
    return srm_out
    
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
                     'Target Gene Name',
                     'Species',
                     'Comments',
                     'Is this a human pluripotent stem cell (hESC or hiPSC) project?'
        ]]
       
       results = srm_df.to_numpy().tolist()
       
       print(srm_df.columns)

       return results

def clicked(a):
    print("clicked") 

def radio_select(self):
    
    radio_var = self.radio_choce
    
    #print(radio_var)
    
    return radio_var

def parse_signature():
    
    try:
        sig_path = os.path.join((os.environ['USERPROFILE']), 'AppData\Roaming\Microsoft\Signatures')
        
        print(sig_path)
        
        os.chdir(sig_path)
        
        
        for file in os.listdir(sig_path):
            #if file.endswith("@stjude.org).htm"):
            if file.startswith('Email signature') and file.endswith(".htm"):
                #print(f"The signature file: {file}")
                sig_htm = file
            else:
                print("No sig found")   
        
        with open(sig_htm, "r") as f:
            html = f.read()
            print(html)
        #match everything between and including the two body tags
        body_pattern = "(?:<style)(.|\n)*?<\/html>"

        #returns match object.  matched text is accessed by .group() because reasons?  Strip newline to condense signature
        text = re.search(body_pattern, html)
        
        print(text)
        
        sig = text.group()
    except:
        sig =""
        
    print(sig)
    return sig
    

if __name__ == "__main__":
  
    xls = 'Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_boostrap_emailer\CAGEServices_Excel Export.xls'   
    #df_from_template(xls)
    
    parse_signature()