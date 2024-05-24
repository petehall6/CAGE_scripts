import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (open_file,
                               df_from_ngs_template,
                               parse_signature
                            )
import pandas as pd
import numpy as np
import shutil
import os
import win32com.client
import glob
import datetime


class NGS_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.button_container = tbs.Frame(self)
        self.button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.user_comments = tbs.StringVar(value="")
        self.ngs_date = tbs.StringVar(value="")

        
        
        self.data = []

        self.create_labels()
        
        self.create_srm_load_btn()
        self.create_ngs_date_field()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        
        self.table = self.create_table()

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "NGS Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
        
        
        self.excel_lbl = tbs.Label(
            master = self.button_container, 
            text="SRM template", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.ngs_lbl = tbs.Label(
            master = self.button_container, 
            text="NGS Date:", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10)
        self.ngs_lbl.grid(column=0, row=2)
    
    def create_ngs_date_field(self):
        
        self.ngs_date_box = tbs.Entry(
        master = self.button_container,
        bootstyle = "info",
        width=22
        )
            
        self.ngs_date_box.grid(column=1, row=2)
        
    def create_srm_load_btn(self):
        
        self.srm_load_btn = tbs.Button(
            master = self.button_container,
            text = "Select SRM Template",
            command = self.load_srm,
            bootstyle=SUCCESS,
            width=25
        )
        
        self.srm_load_btn.grid(column=0,row=1, pady=10)
        

    def create_gen_emails_btn(self):

        self.gen_emails_btn = tbs.Button(
            master = self.button_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.gen_emails_btn.grid(column=0, row=3, pady=10)

    def create_clear_btn(self):

        self.clear_btn = tbs.Button(
            master = self.button_container,
            text = "Clear Entries",
            command = self.clear_controls,
            bootstyle = DANGER,
            width = 25   
        )
        
        self.clear_btn.grid(column=0, row=4,pady=60)
        
        return

    def create_table(self):
        columns = [
            {"text":"SRM Order#"},
            {"text":'Requested By'},
            #{"text":'PI'}, waiting on RIS to add to excel export
            {"text":'Project Number'},
            {"text": "Gene"},
            {"text":'User Comments'},
        ]

        table = Tableview(
            master = self,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=False,
            bootstyle=PRIMARY,
            stripecolor=LIGHT,
        )

        #table.view.selection_set(0)
        #table.view.bind("<<TreeviewSelect>>", clicked(self.table.view.yview()))
        
        table.pack(side=BOTTOM,fill=BOTH, expand=YES, padx=10, pady=10)
        
        
        return table

    def clear_controls(self):
        
        table = self.table
        
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        
        self.data = []

    def load_srm(self):
        self.table.unload_table_data()
        self.data=[]
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        srm_list = df_from_ngs_template(template)
        
        #append table to bottom of frame
        #loop through index
        #then loop through element

        for srm in srm_list:
            for entry in srm:
                self.srm_order = srm[0]
               # self.PI = srm[1], change indices
                self.project_number = srm[1]
                self.requested_by = srm[2]
                self.gene = srm[3]
                self.user_comments = srm[4]

                
            self.data.append((self.srm_order,
                              #self.PI,
                              self.requested_by,
                              self.project_number,
                              self.gene,
                              self.user_comments,
            ))
            
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        def _get_ngs_date(self):
            
            ngs_date = self.ngs_date_box.get()
            
            print(f"NGS Date: {ngs_date}")
            return ngs_date
            
        
        def _get_subject_line(scope, gene, cell_line, objective):
            
            sub_line = "subject"
                    
            return sub_line
        
        #pass email object, cage num
        def _get_attachment(email, srm_order_num):
            #find powerpoint
            try:
                ngs_dir = "Z:\ResearchHome\Groups\millergrp\home\common\\NGS"
                for name in glob.glob(os.path.join(ngs_dir, "*{}".format(srm_order_num))):
                    folder = name

                os.chdir(folder)
                ppt_list = glob.glob("*.pptx")
                latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)
    
            except:
                print("couldn't find slidedeck in CORE Project folder")
                print("Project Number = {}".format(project_num))
                latest_ppt = None

            if latest_ppt is not None:
                email.Attachments.Add(latest_ppt)
                
            
            
            return email
        
        def _body_builder(requester, pi, scope, cell_line, objective, line_lead):
            
            pi = pi.split(", ")[1]
            requester = requester.split(", ")[1]
            
            body=f"""
            <font face="Calibri, Calibri, monospace">
            Hi {pi} and {requester},
            <br><br>
            Great news! Your {gene} {cell_line} NGS data is attached.
            <br><br>
            Best,
            <br><br>
            </font>               
            """
                
            return body
        
        #gets only rows shown in table and acceses those to generate emails
        intact_rows = self.table.get_rows(visible=True)        
        srm_entries=[]
        
        for row in intact_rows:
            srm_entries.append(row.values)
                
        '''
        'SRM Order #',0
        'PI',1 #052024-Waiting on RIS to add PI column.  all other indicies -= 1
        'Requested By',2
        'Project Number',3
        'Gene', 4
        'User Comments',5
        '''    
        sig = parse_signature()
        
        ngs_date = _get_ngs_date(self)

        #self data is a list of list.  loop through each entry to access each field
        for entry in srm_entries:
            
            srm_order_num, requester, project_num, user_comment, gene = entry
            
            #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            email_recip = [requester]
            email_cc = [pi,line_lead]
            
            email_sub = _get_subject_line(scope,gene,pi, )

            email = _get_attachment(email,srm_order_num)

            body = _body_builder(requester,pi,scope,cell_line,objective, line_lead)

            email.To = ";".join(email_recip) #Requester
            email.CC = ";".join(email_cc).replace(".","") #shondra and PI

            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)


if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()