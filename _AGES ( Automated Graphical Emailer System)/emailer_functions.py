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
_PETE = 'Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_boostrap_emailer'

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
                     'Cell Line Lead',
                     'Is this a human pluripotent stem cell (hESC or hiPSC) project?'
        ]]
       
       results = srm_df.to_numpy().tolist()
       #print(results)
       return results

def df_from_design_template(template):
       info = pd.read_excel(template)
       
       tmp = pd.DataFrame(info)

       srm_df = tmp[['SRM Project #',
                     'Project Number',
                     'PI',
                     'Requested By',
                     'Scope',
                     'Target Gene Name',
                     '# of gRNAs?',
                     'gRNA 1',
                     'gRNA 2',
                     'gRNA 3',
                     'gRNA Design - Final PowerPoint File',
        ]]
       
       results = srm_df.to_numpy().tolist()
       #print(results)
       return results

def df_from_tails_template(template):
    
    info = pd.read_excel(template)
    
    tmp = pd.DataFrame(info)
    
    srm_df = tmp[['Principal Investigator',
                  'Entered By',
                  'Date Ordered',
                  'CAGE Project #',
                  'Gene Name/Gene ID',
                  'Number of Tube Samples/Plates',
                  'Sample Format',
                  'Sample Type',
                  'Specify Sample Type',
                  'SRM Order #',
                  'SRM Sample #',
                  'Requested By',
                  'Consolidation Plate?',
                  'Number of Consolidation Plates',
                  'User Comments',
                  'Lab Comments'
                  ]]
    
    
    results = srm_df.to_numpy().tolist()
    
    return results

def df_from_ngs_template(template):
       info = pd.read_excel(template)
       
       tmp = pd.DataFrame(info)

       srm_df = tmp[['SRM Order #',
                     'CAGE Project #',
                     #'PI',
                     'Requested By',
                     'Gene Name/Gene ID',
                     'User Comments'

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

#joke
def window_size():
    user = os.environ.get('USERNAME')
    if user == 'jklein1':
        app = tbs.Window(
            title="CAGE Emailer",
            themename = "vapor",
            size=(1600,800),
            resizable=(True,True),
        )

    else: 
        app = tbs.Window(
            title="CAGE Emailer",
            themename = "superhero",
            size=(1600,800),
            resizable=(True,True),
        )
    print(user)
    return app

if __name__ == "__main__":
  
    xls = 'Z:\ResearchHome\Groups\millergrp\home\common\Python\_pete\_boostrap_emailer\CAGEServices_Excel Export.xls'   
    #df_from_template(xls)
    
    parse_signature()