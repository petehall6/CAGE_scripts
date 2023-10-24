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
        self.button_container = tbs.Frame(self)
        self.button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.project_scope = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.project_objective = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.line_lead = tbs.StringVar(value="")
        
        
        self.data = []

        self.create_labels()
        
        self.create_srm_load_btn()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        
        self.table = self.create_table()

        #TODO FIX the table selection here.  table.view == ttk.Treeview
        #self.table.view.selection
        #self.table.view.bind("<<TreeviewSelect>>", clicked(self.table.view.yview()))
    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Billing Emailer",
            bootstyle = WARNING,
        )
        
        
        self.excel_lbl = tbs.Label(
            master = self.button_container, 
            text="SRM template", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10)
    

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
        
        self.gen_emails_btn.grid(column=0, row=2, pady=10)

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
            {"text":'PI'},
            {"text":'Requested By'},
            {"text":'Project Number'},
            {"text":'Project Scope'},
            {"text":'Cell Line of Choice'},
            {"text":'Project Objective'},
            {"text":'Target Gene Name'},
            {"text":'Cell Line Lead'}
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
        srm_list = df_from_template(template)
        
        #append table to bottom of frame
        #loop through index
        #then loop through element

        try:
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
                    self.species = srm[8]
                    self.line_lead = srm[9].split("Lead-")[1]
                    self.stem_cell = srm[10]
                    
                self.data.append((self.srm_order,
                                self.PI,
                                self.requested_by,
                                self.project_number,
                                self.project_scope,
                                self.cell_line,
                                self.project_objective,
                                self.gene,
                                self.line_lead
                ))
        except: #will catch if no lead is in comments
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
                    self.species = srm[8]
                    self.line_lead = ""
                    self.stem_cell = srm[10]
                    
                self.data.append((self.srm_order,
                                self.PI,
                                self.requested_by,
                                self.project_number,
                                self.project_scope,
                                self.cell_line,
                                self.project_objective,
                                self.gene,
                                self.line_lead
                ))    
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
                latest_ppt = None

            if latest_ppt is not None:
                email.Attachments.Add(latest_ppt)
                
            
            
            return email
        
        def _body_builder(greeting, scope, cell_line, objective, line_lead):
            
            if scope.lower() == "edited cell pool":

                body=f"""{greeting},
                <br><br>
                Great news! Your {gene} {cell_line} edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.
                <br><br>
                The last slide is the most informative.  We were able to get over <font color=red>XX%</font> total editing in the pool with <font color=red>~XX%</font> out of frame indels.
                <br><br>
                We have a contactless pickup system in place.  Please coordinate with {line_lead.split(" ")[0]} to determine a good time window for you to pick up these cells. 
                During the agreed upon time, {line_lead.split(" ")[0]} will place the frozen vials of cells in a dry ice bucket in M4170. 
                The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  
                If you would like the live cultures as well, please come in the next day or so.  
                The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.
                <br><br>
                Don't hesitate to contact me if you have any questions.
                <br><br>
                Best,
                <br><br>
                SM
                <br><br>                
                """
                
            elif scope.lower() == "cell fitness/dependency assay":
                body=f"""{greeting},
                <br><br>
                Great news! Your {gene} {cell_line} fitness assay is complete. Please see the attached slide deck for details.
                <br><br>
                We do/do not see a strong dependency for this gene in this background.
                <br><br>
                Please let me know if you have any questions.
                <br><br>
                Best,
                <br><br>
                SM
                <br><br>
                """
                
            else: 
                body=f"""{greeting},
                <br><br>
                Great news! Your {cell_line} {gene} {objective} project is complete and ready for pick up.  Please see the attached slide deck for details.
                <br><br>
                We have a contactless pickup system in place.  Please coordinate with {line_lead.split(" ")[0]} to determine a good time window for you to pick up these cells. 
                During the agreed upon time, {line_lead.split(" ")[0]} will place the frozen vials of cells in a dry ice bucket in M4170. 
                The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  
                The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.
                <br><br>
                Best,
                <br><br>
                SM
                <br><br>
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
        sig = parse_signature()
        #self data is a list of list.  loop through each entry to access each field
        for entry in srm_entries:
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = entry
            
            #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = [requester,pi]
            email_recip = list(set(recip_list))
            
            if len(email_recip) > 1:
                greeting = f"Hi {pi.split(',')[1]} and {requester.split(',')[1]}"
            else:
                greeting = f"Hi {pi.split(',')[1]}"
            
            email_cc = [line_lead]
                            
            email_sub = _get_subject_line(scope,gene,cell_line, objective)

            email = _get_attachment(email,project_num)

            body = _body_builder(greeting,scope,cell_line,objective, line_lead)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc).replace(".","")

            email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)


if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()