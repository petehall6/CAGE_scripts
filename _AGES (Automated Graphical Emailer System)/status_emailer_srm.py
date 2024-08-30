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

class Status_Tab_srm(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.button_container = tbs.Frame(self)
        self.button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.pi = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.project_scope = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.project_objective = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.line_lead = tbs.StringVar(value="")
        self.weeks = tbs.StringVar(value="")
        self.status_choice = tbs.StringVar(value="")
        self.data = []
        
        self.create_labels()
        self.create_srm_load_btn()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        self.create_textbox()
        self.create_radiobtns()
        self.table = self.create_table()

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Project Status Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
            
        self.excel_lbl = tbs.Label(
            master = self.button_container, 
            text="SRM Template", 
            font=(10), 
            bootstyle = SUCCESS
        )
        
        self.weeks_lbl = tbs.Label(
            master = self.button_container,
            text = 'Number of Weeks: ',
            font = (10),
            bootstyle = 'Success'
        )
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10, columnspan=2)
        self.weeks_lbl.grid(column=0, row=3, pady=10, sticky=W)

    def create_srm_load_btn(self):
        
        self.srm_load_btn = tbs.Button(
            master = self.button_container,
            text = "Select SRM Template",
            command = self.load_srm,
            bootstyle=SUCCESS,
            width=25
        )
        
        self.srm_load_btn.grid(column=0,row=1,pady=10)
        
    def create_radiobtns(self):
        
        self.pool_radiobtn = tbs.Radiobutton(
            master = self.button_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Confirmed Pool",
            value = "Confirmed Pool",
        )
        
        self.screen_radiobtn = tbs.Radiobutton(
            master = self.button_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Initial Screen",
            value = "Initial Screen",
        )
        
        self.delay_screen_radiobtn = tbs.Radiobutton(
            master = self.button_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Delayed Initial Screen",
            value = "Delayed Initial Screen",
        )
        
        self.delay_clone_radiobtn = tbs.Radiobutton(
            master = self.button_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Delayed Clone Hand-off",
            value = "Delayed Clone Hand-off",
        )
        
        self.pool_radiobtn.grid(column=0,row=2,sticky=W, pady=10)
        self.screen_radiobtn.grid(column=1,row=2,sticky=W+E, pady=10)
        self.delay_screen_radiobtn.grid(column=2,row=2,sticky=W, pady=10)
        self.delay_clone_radiobtn.grid(column=3,row=2,sticky=E, pady=10, padx=10)

    def create_gen_emails_btn(self):

        self.gen_emails_btn = tbs.Button(
            master = self.button_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.gen_emails_btn.grid(column=0, row=4, pady=10)

    def create_clear_btn(self):

        self.clear_btn = tbs.Button(
            master = self.button_container,
            text = "Clear Entries",
            command = self.clear_controls,
            bootstyle = DANGER,
            width = 25   
        )
        
        self.clear_btn.grid(column=0, row=5,pady=60)
        
    def create_textbox(self):
        
        self.weeks_box = tbs.Entry(
        master = self.button_container,
        bootstyle = PRIMARY,
        textvariable=self.weeks,
        )
        
        self.srm_label = tbs.Entry(
            
        )
        self.weeks_box.grid(column=1, row=3, sticky='w')
        
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
        
        self.weeks_box.delete(0, 'end')
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
                    self.pi = srm[1]
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
                                self.pi,
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
                    self.pi = srm[1]
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
                                self.pi,
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
        
        #make a way to loop through each project in the srm if multiple lines detected
        
        entries = [self.pi, self.requested_by, self.status_choice.get(), self.weeks.get(), self.gene, self.cell_line, self.project_objective]

        pi, requested_by, status, weeks, gene, cell_line, objective = entries
        
        signature = parse_signature()
        
        def _get_subject_line(gene, cell_line, objective):
        
            sub_line = f"{gene} {objective} {cell_line} status update"
            
            return sub_line

        def _body_builder(greeting, status, weeks, cell_line, objective):

            if status.upper() == "CONFIRMED POOL":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great News! We have successfully confirmed the desired edit in the cell pool for your {cell_line} {gene} {objective} project.  
                We have already sorted for single cells into 96-well plates and will update you when we have screened 
                the plates and identified correctly edited clones. Each modification and cell line is a custom project, 
                and the time will differ widely for each project depending on several factors. Based on the details of your 
                specific project, we estimate that we will have identified clones in about 12 weeks. 
                If you have any questions or concerns, please don't hesitate to reach out.
                <br><br>
                </font>
                """
                
            elif status.upper() == "INITIAL SCREEN":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great News! We have successfully identified clones with the desired modification for your {cell_line} {gene} {objective} project. 
                If everything goes as planned, we expect to hand off these clones to you in 4 weeks. 
                Please let me know if you have any questions.
                <br><br>
                </font>
                """

            elif status.upper() == "DELAYED INITIAL SCREEN": 
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                I wanted to provide you with an update on your {cell_line} {gene} {objective} project. Unfortunately, we were unable to identify any correctly 
                edited clones during the initial screen.  We are reviewing our data and reevaluating the editing strategy now.  
                We are still working hard to obtain the desired edited clone(s), but there is going to be a delay in the timeline 
                as we restart the process.  Please let me know if you have any questions.
                <br><br>
                </font>
                """
                
            elif status.upper() == "DELAYED CLONE HAND-OFF": 
                body=f"""
                <font face="Calibri, Calibri, monospace">{greeting},
                <br><br>
                I wanted to provide you with an update on your {cell_line} {gene} {objective} project.  The clones are growing slower than expected, 
                which has resulted in a delayed timeline for project completion.  We will email you when they are expanded, fully QC’d,
                and ready for pick up. Based on their current rate of growth, I would expect them to be ready in {weeks} weeks.  
                Please let me know if you have any questions.
                <br><br>
                </font>
                """

            return body
        
        #mail generator
        outlook = win32com.client.Dispatch("Outlook.Application")
        email = outlook.CreateItem(0)
        
        #removes duplicates and rephrases the greeting to a single person
        recip_list = [requested_by,pi]
        email_recip = list(set(recip_list))
        
        if len(email_recip) > 1:
            greeting = f"Hi {pi.split(',')[1]} and {requested_by.split(',')[1]}"
        else:
            greeting = f"Hi {pi.split(',')[1]}"
        
        email_cc = "Miller, Shondra"
                        
        email_sub = _get_subject_line(gene,cell_line, objective)

        body = _body_builder(greeting, status, weeks, cell_line, objective)

        email.To = ";".join(email_recip)
        email.CC = email_cc
        if status == "Initial Screen" or status == "Delayed Clone Hand-off":
            email.BCC = "Porter, Shaina"

        email.Subject = email_sub

        #find html signature file in each individual userprofile
        email.HTMLBody = body + signature
        #Display(False) loads all emails at once and gives focus back to ttk window
        email.Display(False)

        return

if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()