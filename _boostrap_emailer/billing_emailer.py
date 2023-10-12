import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (open_file,
                               df_from_template,
                               clicked,
                               parse_signature
                            )
import pandas as pd
import numpy as np
import shutil
import os
import win32com.client


class Billing_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        #self.colors = master_window.style.colors
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.project_scope = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.project_objective = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        
        
        self.data = []
        
        self.create_srm_load_btn()
        self.create_labels()

        
        self.table = self.create_table()
        self.create_gen_emails_btn()
        self.create_clear_btn()

        #TODO FIX the table selection here.  table.view == ttk.Treeview
        #self.table.view.selection
        #self.table.view.bind("<<TreeviewSelect>>", clicked(self.table.view.yview()))
        


    def create_srm_load_btn(self):
        button_container = tbs.Frame(self)
        button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.srm_load_btn = tbs.Button(
            master = button_container,
            text = "Select SRM Template",
            command = self.load_srm,
            bootstyle=SUCCESS,
            width=25
        )
        
        self.srm_load_btn.pack(side=LEFT, padx=5)

    def create_gen_emails_btn(self):
        button_container = tbs.Frame(self)
        button_container.pack(fill=X, expand=YES, pady=(20,10))
        
        self.gen_emails_btn = tbs.Button(
            master = button_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.gen_emails_btn.pack(side=RIGHT, padx=5)

    def create_clear_btn(self):
        button_container = tbs.Frame(self)
        button_container.pack(fill=X, expand=YES, pady=(30))
        
        self.clear_btn = tbs.Button(
            master = button_container,
            text = "Clear Entries",
            command = self.clear_controls,
            bootstyle = DANGER,
            width = 25   
        )
        
        self.clear_btn.pack(side=RIGHT, padx=5)
        
        
        
        return

    def generate_emails(self):
        
        #USERS/**/AppData/Roaming/Microsoft/Signatures
        
        #sig_path = os.path.join((os.environ['USERPROFILE']), 'APPDATAAppData\Roaming\Microsoft\Signatures')        
        
        #not great practice to repack self.data but makes me life easier right now
        
        srm_entries = self.data
                
        for entry in srm_entries:
            
            email_recip = str(entry[2])
            print(f"Primary recipient: {email_recip}")
            email_cc = str(entry[1])
            print(f"CC list: {email_cc}")
            
            body = (
            f"""
            Hi,
            
            The CAGE has done the things we did that you wanted us to do.
            <br>
            <br>
            """
            
            )

            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            #splits each element from list and seperates each with ';' allowing for multiple recipenants and cc's
            email.To = email_recip
            email.CC = email_cc
            email.Subject = "BILLING EMAILER"
            #explicitly converting slashes for ease
            #attachment = os.path.join(DEST_DIR,cage_num+"_fastq_zipped.zip").replace("\\","/")
            #email.Attachments.Add(attachment)
            
            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            email.Display(True)

    def create_labels(self):
        lbl_container = tbs.Frame(self)
        lbl_container.pack(fill=X, expand=YES, pady=5)
        self.excel_lbl = tbs.Label(lbl_container, text="Excel Name", font=(10), bootstyle = "SUCCESS")
        
        self.excel_lbl.pack(side=LEFT, padx=5)

    def create_table(self):
        columns = [
            {"text":"SRM Order#"},
            {"text":'PI'},
            {"text":'Requested By'},
            {"text":'Project Number'},
            {"text":'Project Scope'},
            {"text":'Cell Line of Choice'},
            {"text":'Project Objective'},
            {"text":'Target Gene Name'}
        ]

        table = Tableview(
            master = self,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=True,
            bootstyle=PRIMARY,
            stripecolor=LIGHT,
        )

        #table.view.selection_set(0)
        #table.view.bind("<<TreeviewSelect>>", clicked(self.table.view.yview()))
        
        
        table.pack(side=BOTTOM,fill=BOTH, expand=YES, padx=10, pady=10)
        
        
        return table

    def load_srm(self):
        
        self.data=[]
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        srm_list = df_from_template(template)
        
        print(f"This is the srm_order column: {srm_list}")
        
        print(self.data)
        #append table to bottom of frame
                
        #loop through index
        #then loop through element
        for srm in srm_list:
            for entry in srm:
                self.srm_order = srm[0]
                self.PI = srm[1]
                self.requested_by = srm[2]
                self.project_number = srm[3]
                self.project_scope = srm[4]
                self.cell_line = srm[5]
                self.project_objective = srm[6]
                self.gene = srm[7]
                
            self.data.append((self.srm_order,self.PI,self.requested_by,self.project_number,self.project_scope,self.cell_line,self.project_objective,self.gene))
            
            
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    
        
    def clear_controls(self):
        
        table = self.table
        
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        
        self.data = []
        
        print(self.data)
        
if __name__ == "__main__":
    import emailer_gui
    emailer_gui.app.mainloop()