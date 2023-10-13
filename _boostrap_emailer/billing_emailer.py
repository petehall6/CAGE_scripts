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
import glob
import datetime


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

    def clear_controls(self):
        
        table = self.table
        
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        
        self.data = []

    def load_srm(self):
        
        self.data=[]
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        srm_list = df_from_template(template)
        
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
            
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        def _get_subject_line(scope, gene, cell_line, objective):
            
            if scope.upper() == "EDITED CELL POOL":
                sub_line = f"{gene} {cell_line} Edited Cell Pool Complete"
            elif scope.upper() == "CELL LINE CREATION":
                sub_line = f"{cell_line} {gene} {objective} Cell Line Complete"
            elif scope.upper() == "CELL FITNESS/DEPENDENCY ASSAY":
                sub_line = f"CelFi Assay for {gene} in {cell_line} Cells Complete"
                    
            return sub_line
        
        #pass email object, cage num
        def _get_attachment(email, project_num):
            #find powerpoint
            try:
                path = "Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS/"
                for name in glob.glob(os.path.join(path, "*{}".format(project_num))):
                    folder = name

                os.chdir(folder)
                ppt_list = glob.glob("*.pptx")
                latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)
    
            except:
                print("couldn't find slidedeck in CORE Project folder")
                print("Project Number = {}".format(project_num))

            if latest_ppt is not None:
                email.Attachments.Add(latest_ppt)
            
            
            return email
        
        def _body_builder(requester, pi, scope, cell_line, objective):
            
            pi = pi.split(", ")[1]
            requester = requester.split(", ")[1]
            
            
            if scope.lower() == "edited cell pool":

                body=f"""Hi {pi} and {requester},
                <br><br>
                Great news! Your {gene} {cell_line} edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.
                <br><br>
                The last slide is the most informative.  We were able to get over <font color=red>XX%</font> total editing in the pool with <font color=red>~XX%</font> out of frame indels.
                <br><br>
                We have a contactless pickup system in place right now.  Please coordinate with <b><font color=red>XXXXX</font></b> to let <b><font color=red>XXXX</font></b> know a good time window for you to pick up these cells. 
                During the agreed upon time, <b><font color=red>XXXX</font></b> will place the frozen vials of cells in a dry ice bucket in M4170. 
                The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  
                If you would like the live cultures as well, please come in the next day or so.  
                The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.
                <br><br>
                Don't hesitate to contact me if you have any questions.
                <br><br>
                Best,
                <br><br>
                SM                
                """
                
            elif scope.lower() == "cell fitness/dependency assay":
                body=f"""Hi {pi} and {requester},
                <br><br>
                Great news! Your {gene} {cell_line} fitness assay is complete. Please see the attached slide deck for details.
                <br><br>
                We <font color=red>do/do</font> not see a strong dependency for this gene in this background.
                <br><br>
                Please let me know if you have any questions.
                <br><br>
                Best,
                <br><br>
                SM
                """
                
            else: 
                body=f"""Hi {pi} and {requester},
                <br><br>
                Great news! Your {cell_line} {gene} {objective} project is complete and ready for pick up.  Please see the attached slide deck for details.
                <br><br>
                We currently have a contactless pickup system in place.  Please arrange a time window with <font color=red>XXXXX</font> in which someone can pick up the cells.  
                At the agreed upon time, <font color=red>he/she</font> will place your frozen vials of cells into a dry ice bucket in M4170.  
                The dry ice bucket will be on the counter in front of you as you walk in.  
                Your live cultures will be in the first incubator to the right (top incubator, bottom shelf) and labeled accordingly. Please also bring dry ice for the pickup.
                <br><br>
                As always, please let me know if you have any questions.
                <br><br>
                Best,
                <br><br>
                SM
                """
                
            return body
        
        #gets only rows shown in table and acceses those to generate emails
        intact_rows = self.table.get_rows(visible=True)        
        srm_entries=[]
        
        for row in intact_rows:
            srm_entries.append(row.values)
                
        '''
        'SRM Order #',0
        'PI',1
        'Requested By',2
        'Project Number',3
        'Project Scope',4
        'Cell Line of Choice',5
        'Project Objective',6
        'Target Gene Name' 7
        '''    
        
        #self data is a list of list.  loop through each entry to access each field
        for entry in srm_entries:
            
  
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene = entry
            
            #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            email_recip = str(requester)
            email_cc = str(pi)
            
            email_sub = _get_subject_line(scope,gene,cell_line, objective)

            email = _get_attachment(email,project_num)

            body = _body_builder(requester,pi,scope,cell_line,objective)

            email.To = email_recip
            email.CC = email_cc

            email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)


if __name__ == "__main__":
    import emailer_gui
    emailer_gui.app.mainloop()