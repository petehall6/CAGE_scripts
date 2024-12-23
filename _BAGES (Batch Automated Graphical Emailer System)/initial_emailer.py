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


class Initial_Tab(tbs.Frame):
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

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Initial Batch Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
        
        
        self.excel_lbl = tbs.Label(
            master = self.button_container, 
            text="SRM Template", 
            font=(10), 
            bootstyle = SUCCESS
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

    def create_table(self):
        columns = [
            {"text":'Project Number'},
            {"text":'PI'},
            {"text":'Requested By'},
            {"text":'Project Scope'},
            {"text":'Species'},
            {"text":'Project Objective'},          
            {"text":'Gene'},
            {"text":'gRNA 1'},
            {"text":'gRNA 2'},
            {"text":'gRNA 3'},
            {"text":'gRNA 4'},
            {"text":'PD name'},
            {"text":'ssODN#1 name'},
            {"text":'ssODN#2 name'}
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
        for srm in srm_list:
            for entry in srm:
                self.project_number = srm[0]
                self.PI = srm[1]
                self.requested_by = srm[2]
                self.project_scope = srm[3]
                self.species = srm[4]
                self.project_objective = srm[5]
                self.gene = srm[6]
                self.grna_1 = srm[7]
                self.grna_2 = srm[8]
                self.grna_3 = srm[9]
                self.grna_4 = srm[10]
                self.pd_name = srm[11]
                self.ssodn1_name = srm[12]
                self.ssodn2_name = srm[13]
                
            self.data.append((self.project_number,
                            self.PI,
                            self.requested_by,
                            self.project_scope,
                            self.species,
                            self.project_objective,
                            self.gene,
                            self.grna_1,
                            self.grna_2,
                            self.grna_3,
                            self.grna_4,
                            self.pd_name,
                            self.ssodn1_name,
                            self.ssodn2_name
            ))
            
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        signature = parse_signature()
        
        def _email_writer_single(project_details):
            
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = project_details
            
            def _get_subject_line(scope, gene, cell_line, objective):
            
                if scope.upper() == "EDITED CELL POOL":
                    sub_line = f"{gene} {cell_line} Edited Cell Pool Complete"
                elif scope.upper() == "CELL LINE CREATION":
                    sub_line = f"{cell_line} {gene} {objective} Cell Line Complete"
                elif scope.upper() == "CELL FITNESS/DEPENDENCY ASSAY":
                    sub_line = f"CelFi Assay for {gene} in {cell_line} Cells Complete"    
                return sub_line
            
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

                    body=f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Great news! Your {gene} {cell_line} edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.
                    <br><br>
                    The last slide is the most informative.  We were able to get over <font color=red>XX%</font> total editing in the pool with <font color=red>~XX%</font> out of frame indels.
                    <br><br>
                    We have a contactless pickup system in place.  Please coordinate with {line_lead.split(" ")[0]} to determine a good time window for you to pick up these cells. 
                    During the agreed upon time, {line_lead.split(" ")[0]} will place the frozen vials of cells in a dry ice bucket in M4170. 
                    The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  
                    If you would like the live cultures as well, please come in the next day or so.  
                    Your live cultures will be on the bottom shelf of the "Pick-up" incubator, which is labeled accordingly.  Please bring dry ice for the pickup.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>
                    </font>                
                    """
                    
                elif scope.lower() == "cell fitness/dependency assay":
                    body=f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Great news! Your {gene} {cell_line} fitness assay is complete. Please see the attached slide deck for details.
                    <br><br>
                    We do/do not see a strong dependency for this gene in this background.
                    <br><br>
                    Please let me know if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>
                    </font>
                    """
                    
                else: 
                    body=f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Great news! Your {cell_line} {gene} {objective} project is complete and ready for pick up.  Please see the attached slide deck for details.
                    <br><br>
                    We have a contactless pickup system in place.  Please coordinate with {line_lead.split(" ")[0]} to determine a good time window for you to pick up these cells. 
                    During the agreed upon time, {line_lead.split(" ")[0]} will place the frozen vials of cells in a dry ice bucket in M4170. 
                    The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  
                    Your live cultures will be on the bottom shelf of the "Pick-up" incubator, which is labeled accordingly.  Please bring dry ice for the pickup.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>
                    </font>
                    
                    """
                    
                return body
            
            
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
            
            email.HTMLBody = body + signature
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)

        def _email_writer_multi(project_df):
            
            projects = project_df.values.tolist()
            
            #initizlie all the list at once
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = ([] for i in range(9))
            #unpacked nested list into individual list
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = map(list,zip(*projects))

            #remove duplicates form line_lead list
            line_lead = list(set(line_lead))
            #if more than one project lead, set line_lead to 'XXXXXXX' so Shondra can designate
            if len(line_lead) > 1:
                line_lead.clear()
                line_lead.append("XXXXXXXXX")
            
            def _get_attachment(email, project_num):
                #find powerpoint
                
                for proj in project_num:
                
                    try:
                        path = "Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS/"
                        for name in glob.glob(os.path.join(path, "*{}".format(proj))):
                            folder = name

                        os.chdir(folder)
                        ppt_list = glob.glob("*.pptx")
                        latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)
            
                    except:
                        print("couldn't find slidedeck in CORE Project folder")
                        print("Project Number = {}".format(proj))
                        latest_ppt = None

                    if latest_ppt is not None:
                        email.Attachments.Add(latest_ppt)
                        
                return email
        
            def _bullet_maker(srm_order_num, gene, scope, objective, cell_line):
                bullet_list =""
                for order, proj_gene, proj_scope, proj_obj,cells in zip(srm_order_num,gene,scope, objective, cell_line):
                    print(f"order {order}, proj_gene: {proj_gene}, proj_scope: {proj_scope}, proj_obj: {proj_obj}, cells: {cells}")
                    
                    
                    bullet_list += (f"<li>SRM: {order}- {cells} {proj_gene} {str(proj_obj).replace('Gene','').replace('nan','')} {proj_scope} </li>")
                
                print(f"The bullet_list {bullet_list}")
                    
                return bullet_list
            
            def _body_builder_multi(srm_order_num,greeting, scope, cell_line, objective, line_lead):
                
                if 'XXXX' not in line_lead[0]:
                    line_lead = line_lead[0].split(" ")
                
                
                bullets = _bullet_maker(srm_order_num,gene,scope, objective, cell_line)
                
                body = f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news! The following projects are ready for pickup.  Please see the attached slide decks for details:
                <br><br>
                <ul>
                {bullets}
                </ul>
                <br>
                We have a contactless pickup system in place. Please arrange a time window with {line_lead[0]} in which someone can pick up the cells.  
                At the agreed upon time, {line_lead[0]} will place your frozen vials of cells into a dry ice bucket in M4170. The dry ice bucket will 
                be straight in front as you walk in. Your live cultures will be on the bottom shelf of the "Pick-up" incubator, which is labeled accordingly.
                <br><br>
                As always, please let me know if you have any questions.<br>
                <br>
                Best,<br>
                SM
                <br><br>
                </font>

                """

                return body
            
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = requester + pi
            
            email_recip = list(set(recip_list))
            
            print(f"Email receip list: {email_recip}")
            
            if len(email_recip) > 1:
                first_names=[]
                for req in requester:
                    first_names.append(req.split(',')[1])

                first_names = list(set(first_names))
                #insert 'and' in front of last element            
                first_names[-1] = ' and '+first_names[-1] 
     
                greeting = f"Hi {pi[0].split(',')[1]}, {(','.join(first_names))}"
                
            else:
                greeting = f"Hi {pi[0].split(',')[1]},"
            
            if 'XXXX' in line_lead:
                email_cc =""
            else:
                email_cc = line_lead
                            
            email_sub = "Projects ready for pickup"

            email = _get_attachment(email,project_num)

            body = _body_builder_multi(srm_order_num,greeting,scope,cell_line,objective, line_lead)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc)

            email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            
            email.HTMLBody = body + signature
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)
            
            
            return
        #gets only rows shown in table and access those to create the emails
        intact_rows = self.table.get_rows(visible=True)
        srm_entries=[]
        for row in intact_rows:
            srm_entries.append(row.values)

        #convert to df to create pi_specific sub frame
        columns = [
                     'SRM Order #',
                     'PI',
                     'Requested By',
                     'Project Number',
                     'Project Scope',
                     'Cell Line of Choice',
                     'Project Objective',
                     'Target Gene Name',
                     'Line Lead',
        ]
        
        srm_entries_df = pd.DataFrame(srm_entries, columns=columns)
        
        print(srm_entries_df)
        
        #get a list of unique PIs in request
        pi_list = list(set(srm_entries_df["PI"].values.tolist()))

        #loop create pi specific df by loop though srm_df matching for pi name

        for pi in pi_list:
            pi_specifc_df = srm_entries_df.loc[srm_entries_df['PI'] == pi]

            #check for number of rows and pass to either multi or single project format
            #multiple entries per PI
            if pi_specifc_df.shape[0] > 1:
                    #print(f"Multiple projects: {pi_specifc_df.iloc[0][1]}")
                    #pass to body_builder_multi
                    _email_writer_multi(pi_specifc_df) 
            else:
                    #print(f"single project: {pi_specifc_df.iloc[0][1]}")
                    #convert back into list and pass to writer_single
                    #kept in list form since single mode was originally written and multi project was an added on feature
                    proj_details = pi_specifc_df.values.tolist()[0]
                    
                    _email_writer_single(proj_details)
                    
        return
        

if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()