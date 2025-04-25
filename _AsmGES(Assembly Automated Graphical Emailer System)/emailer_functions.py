from tkinter import filedialog
from tkinter import ttk
import ttkbootstrap as tbs
import os
import pandas as pd
import numpy as np
import itertools
import codecs
import re

DESKTOP = os.path.join(os.path.join(os.environ['USERPROFILE'], 'Downloads'))

#functions
def open_file():
    
    fileytpes = (
        ('Excel files', '*.xls'),
        ('All files', '*.*')
    )
    
    srm_out = filedialog.askopenfilename(
        title = "Select file",
        #initialdir= _PETE,
        initialdir= DESKTOP,
        filetypes=fileytpes  
    )
    
    return srm_out
    
def df_from_template(template):
    cols = [
        'SRM Order #',
        'Project Number',
        'PI',
        'Requested By',
        'Scope',
        'Gene',
        'Species',
        'Specify Vector Preference',
        'gRNA 1',
        'gRNA 2',
    ]
    info = pd.read_excel(template, usecols=cols)
    
    tmp = pd.DataFrame(info)
        
    srm_df = tmp[[
        'SRM Order #',
        'Project Number',
        'PI',
        'Requested By',
        'Scope',
        'Gene',
        'Species',
        'Specify Vector Preference',
        'gRNA 1',
        'gRNA 2',
    ]]
    
    results = srm_df.to_numpy().tolist()
    #print(results)
    return results


def clicked(a):
    print("clicked") 

def radio_select(self):
    
    radio_var = self.radio_choice
    
    return radio_var

def parse_signature():
    
    try:
        sig_path = os.path.join((os.environ['USERPROFILE']), 'AppData\Roaming\Microsoft\Signatures')
        
        print(sig_path)
        
        os.chdir(sig_path)
        
        for file in os.listdir(sig_path):
            #if multiple .htm signatures it will just overwrite.  Will need to fix on individual basis
            if file.endswith(".htm"):
                sig_htm = file
                
        if sig_htm is not None:
            print(f"Signature found {sig_htm}")

            with codecs.open(sig_htm, "r", encoding='utf-8',
                            errors='ignore') as f:
                html = f.read()
            f.close()
            #match everything between and including the two body tags
            body_pattern = "(?:<style)(.|\n)*?<\/html>"

            #returns match object.  matched text is accessed by .group() because reasons?  Strip newline to condense signature
            text = re.search(body_pattern, html)

            sig = text.group()
            
            open_tag = '<font face="Calibri, Calibri, monospace">'
            close_tag = '</font>'
            
            sig = open_tag+sig+close_tag
            
            
        else:
            sig = " "
    except:
        print("No signature found.  I will still work though...hopefully")
        sig = " "
    return sig


if __name__ == "__main__":
    #xls = 'Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_boostrap_emailer\CAGEServices_Excel Export.xls'   
    #df_from_template(xls)
    
    parse_signature()