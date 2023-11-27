import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (open_file,
                               df_from_tails_template,
                               parse_signature,
                            )
import pandas as pd
import numpy as np
import shutil
import os
import win32com.client
import glob
import datetime


class Tails_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.header_container = tbs.Frame(self)
        self.header_container.pack(side=TOP,fill=X, expand=YES, pady=(15,10))
        
        self.table_container = tbs.Frame(self)
        self.table_container.pack(fill=BOTH,expand=YES, pady=(5,5))
        
        self.button_container = tbs.Frame(self)
        self.button_container.pack(side=BOTTOM, fill=X, expand=YES, pady=(10,10))
        
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
        self.success_choice = tbs.StringVar(value="")
        self.edit_choice = tbs.StringVar(value="")
    
        
        self.data = []

        self.create_labels()
        self.create_buttons()
        self.create_comboxes()

        self.table = self.create_table()

    def create_buttons(self):
        
        self.srm_load_btn = tbs.Button(
            master = self.header_container,
            text = "Select SRM Template",
            command = self.load_srm,
            bootstyle=SUCCESS,
            width=25
        )
        
        self.gen_emails_btn = tbs.Button(
            master = self.button_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.clear_btn = tbs.Button(
            master = self.button_container,
            text = "Clear Entries",
            command = self.clear_controls,
            bootstyle = DANGER,
            width = 25
        )
        
        self.store_btn = tbs.Button(
            master = self.button_container,
            text = 'Store',
            command =self.store_values,
            bootstyle = SUCCESS,
            width = 25
        )
        
        self.next_btn = tbs.Button(
            master = self.button_container,
            text = 'Next',
            command = self.nextbtn_click,
            bootstyle = INFO,
            width = 25
        )
        
        self.prev_btn = tbs.Button(
            master = self.button_container,
            text = 'Previous',
            command =self.prvbtn_click,
            bootstyle = INFO,
            width = 25
        )
        
        self.srm_load_btn.grid(column=0,row=1, pady=10)
        self.store_btn.grid(column=1, row=6,pady=10, padx=5)
        self.prev_btn.grid(column=3, row=6,pady=10,padx=5)
        self.next_btn.grid(column=4, row=6,pady=10,padx=5)
        self.gen_emails_btn.grid(column=0, row=7, pady=60)
        self.clear_btn.grid(column=2, row=7,pady=60, sticky=E, padx = 10)
        
    def create_labels(self):    
        self.title_lbl = tbs.Label(
            master = self.header_container,
            text = "Tails Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
        
        self.excel_lbl = tbs.Label(
            master = self.header_container, 
            text="SRM Template", 
            font=(10), 
            bootstyle = SUCCESS
        )
        
        self.proj_lbl = tbs.Label(
            master = self.button_container,
            text = "Current Project: ",
            font=(10),
            bootstyle = SUCCESS
        )
        
        self.active_proj_lbl = tbs.Label(
            master = self.button_container,
            text = 'PLACE HOLDER',
            font=(10),
            bootstyle = SUCCESS
        )
        
        self.success_lbl = tbs.Label(
            master = self.button_container,
            text = 'Success?',
            font=(10),
            bootstyle = SUCCESS
        )
        
        self.edit_lbl = tbs.Label(
            master = self.button_container,
            text = 'Edit:',
            font=(10),
            bootstyle = SUCCESS
        )
        
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10, sticky=W)
        self.proj_lbl.grid(column = 0,row=3, sticky=E)
        self.active_proj_lbl.grid(column=1,row=3, columnspan=3, sticky=W+E+N+S)
        self.success_lbl.grid(column=0,row=4, pady=10)
        self.edit_lbl.grid(column=0, row=5, pady=10)
    
    def create_table(self):
        columns = [
            {"text":"SRM Order#"},
            {"text":'PI'},
            {"text":'Requested By'},
            {"text":'Entered By'},
            {"text":'Project Number'},
            {"text":'Gene'},
            {"text":'Number of Sample'},
            {"text":'Sample Format'},
            {"text":'Sample Type'},
            {"text":'Success'},
            {"text":'Edit'},
        ]

        self.table = Tableview(
            master = self.table_container,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=False,
            bootstyle=PRIMARY,
            stripecolor=LIGHT,
            autoalign=False,
        )
        
        #self.table.view.selection_set(0)
        self.table.view.bind("<<TreeviewSelect>>", self.tableview_clicked)
        
        self.table.pack(side=BOTTOM,fill=BOTH, expand=YES, padx=10, pady=10)
        
        return self.table

    def create_comboxes(self):
        
        
        outcomes = [' ','Yes','No']
        edits = [' ','KO','KI','CKO','Del','ssODN','PM','Data']
        
        self.success_box = tbs.Combobox(
            master = self.button_container,
            bootstyle = "info",
            value = outcomes,
        )
        
        self.edits_box = tbs.Combobox(
            master = self.button_container,
            bootstyle = "info",
            value = edits,
        )

        self.success_box.grid(column=1,row=4, sticky=W)
        self.edits_box.grid(column=1, row=5, sticky=W)





    def load_srm(self):
        self.table.unload_table_data()
        self.data=[]
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        srm_list = df_from_tails_template(template)
        
        #append table to bottom of frame
        #loop through index
        #then loop through element
            
        for srm in srm_list:
            for entry in srm:
                self.srm_order = srm[9]
                self.PI = srm[0]
                self.requested_by = srm[10]
                self.entered_by = srm[1]
                self.project_number = srm[3].strip().split(" and")[0]
                self.gene = srm[4]
                self.sample_num = srm[5]
                self.sample_format = srm[6]
                self.sample_type = srm[7]                   
            self.data.append((self.srm_order,
                            self.PI,
                            self.requested_by,
                            self.entered_by,
                            self.project_number,
                            self.gene,
                            self.sample_num,
                            self.sample_format,
                            self.sample_type
            ))
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def tableview_clicked(self,event):
        
        #item will return a tuple with text,values and other info.  access info
        selected_proj_info = self.table.view.item(self.table.view.focus(),"values"[0])
        
        #set label to PI and srm#
        pi = selected_proj_info[1].split(",")[0]
        srm_no = selected_proj_info[0]
        self.active_proj_lbl.configure(text=f"SRM: {srm_no}. PI: {pi}")
        
        if len(selected_proj_info) != 9:
            print("setting combo value")
            self.success_box.set(selected_proj_info[9])
            self.edits_box.set(selected_proj_info[10])
        else:
            self.success_box.set(" ")
            self.edits_box.set(" ")

    def generate_emails(self):
        
        signature = parse_signature()
        initial_choice = self.initial_choice.get()
        
        def _email_writer_single(project_details,initial_choice):
            
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = project_details
            
            def _get_subject_line(gene, srm_order_num):
            
                sub_line = f"gRNA Designs for {gene}, SRM order {srm_order_num}"
                
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
        
            def _body_builder(greeting, gene, scope, initial_choice):

                body=f"""{greeting},
                <br><br>
                Great news!  Attached are the designs, off-target analysis, and our recommendations for which gRNAs to move forward 
                with for your {gene} ({scope}).
                <br><br>
                Please let me know if you have any questions or if you would like to move forward with our recommendations.
                <br><br>
                Best,
                <br><br>
                {initial_choice}
                <br><br>
                """
                    
                return body
            
            
        #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = [requester,pi]
            email_recip = list(set(recip_list))
            #copies either Jon or Barnanda depending on who sent the email
            email_cc = ["Shondra Miller"]
            if initial_choice == "JK":
                email_cc.append("Baranda Hanson")
            else:
                email_cc.append("Jon Klein")
            
            
            if len(email_recip) > 1:
                greeting = f"Hi {pi.split(',')[1]} and {requester.split(',')[1]}"
            else:
                greeting = f"Hi {pi.split(',')[1]}"
                            
            email_sub = _get_subject_line(scope,srm_order_num)

            email = _get_attachment(email,project_num)

            body = _body_builder(greeting,gene,scope,initial_choice)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc).replace(".","")
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            email.HTMLBody = body + signature
            
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)

        def _email_writer_multi(project_df, initial_choice):
            
            projects = project_df.values.tolist()
            
            #initizlie all the lsit at once
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = ([] for i in range(9))
            #unpacked nested list into individual list
            srm_order_num, pi, requester, project_num, scope, cell_line, objective, gene, line_lead = map(list,zip(*projects))

            def _get_attachment(email, project_num):
                #find powerpoint
                
                for proj in project_num:
                
                    try:
                        path = "Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS/"
                        for name in glob.glob(os.path.join(path, "*{}".format(proj))):
                            folder = name
                        print(f"folder: {folder}")
                        os.chdir(folder)
                        ppt_list = glob.glob("*.pptx")
                        latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)
            
                    except:
                        print("couldn't find slidedeck in CORE Project folder")
                        print(f"Project Number = {proj}")
                        latest_ppt = None

                    if latest_ppt is not None:
                        email.Attachments.Add(latest_ppt)
                        
                return email
        
            def _bullet_maker(srm_order_num, gene, scope):
                bullet_list =""
                for order, proj_gene, proj_scope in zip(srm_order_num,gene,scope):
                    bullet_list += (f"<li>SRM: {order}: {proj_gene} ({proj_scope})</li>")
                    
                print(f"The bullet_list {bullet_list}")
                    
                return bullet_list

            def _body_builder(srm_order_num,greeting, scope, initial_choice):
                
                bullets = _bullet_maker(srm_order_num,gene,scope)
                
                body = f"""{greeting},
                <br><br>
                Great news! Attached are the designs, off-target analysis, 
                and our recommendations for which gRNAs to move forward with for the following projects:
                <br><br>
                <ul>
                {bullets}
                </ul>
                Please let me know if you have any questions or if you would like to move forward with our recommendations.
                <br><br>
                Best,<br>
                {initial_choice}
                <br><br>
                """

                return body
            
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = requester + pi
            
            email_recip = list(set(recip_list))
            
            if len(email_recip) > 1:
                first_names=[]
                for req in requester:
                    first_names.append(req.split(',')[1])
    
                #insert 'and' in front of last element            
                first_names[-1] = ' and '+first_names[-1] 
     
                greeting = f"Hi {pi[0].split(',')[1]}, {(','.join(first_names))}"
            else:
                greeting = f"Hi {pi[0].split(',')[1]}"
            
             #copies either Jon or Barnanda depending on who sent the email
            email_cc = ["Shondra Miller"]
            if initial_choice == "JK":
                email_cc.append("Baranda Hanson")
            else:
                email_cc.append("Jon Klein")
                            
            email_sub = "Projects ready for pickup"

            email = _get_attachment(email,project_num)

            body = _body_builder(srm_order_num,greeting,scope,initial_choice)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc)

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
                    _email_writer_multi(pi_specifc_df,initial_choice) 
            else:
                    #print(f"single project: {pi_specifc_df.iloc[0][1]}")
                    #convert back into list and pass to writer_single
                    #kept in list form since single mode was originally written and multi project was an added on feature
                    proj_details = pi_specifc_df.values.tolist()[0]
                    
                    _email_writer_single(proj_details,initial_choice)
                    
        return

    def clear_controls(self):
        
        table = self.table
        
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        
        self.success_box.set("")
        self.edits_box.set(" ")
        
        
        self.data = []

    def store_values(self):
        
        #get combobox input
        self.success_choice = self.success_box.get()
        self.edit_choice = self.edits_box.get()

        #get row tuple where focus is (highlighted row)
        selected_proj_info = self.table.view.item(self.table.view.focus(),"values")
        
        proj_appended=[]
        
        #transfer tuple values to list for appending and any necessary overwriting
        #if line has already been edited (len(tuple)=9), replace edit and success indecies
        if len(selected_proj_info) == 9:
            for info in selected_proj_info:
                proj_appended.append(info)
                
        else:#only get first 9 tuple values
            for info in selected_proj_info[:9]:
                proj_appended.append(info)
        
        proj_appended.append(self.success_choice)
        proj_appended.append(self.edit_choice)
        
        #convert back to tuple to plug back into item()
        updated_proj_info = tuple(proj_appended)
        #text is normally blank, have to specifiy value update
        self.table.view.item(self.table.view.focus(),text="",values=updated_proj_info)

    def nextbtn_click(self):
        #get current table index
        curr_index = int(self.table.view.focus()[1:])
        print(curr_index)
        
        #add 1 and any necessary 0's to make it 00n.  
        next_index = str("I"+(str(curr_index+1).zfill(3)))
        print(f"next index: {next_index}")
        #loop back through list
        try:
            self.table.view.selection_set(next_index)
            self.table.view.focus(next_index)
            print(f"Moved focus: {next_index}")
        except:#go back to the top of the list if it overruns num of rows
            self.table.view.selection_set('I001')
            self.table.view.focus('I001')
    
    def prvbtn_click(self):
        #get current table index, strip 'I'
        curr_index = int(self.table.view.focus()[1:])
        
        #add 1 and any necessary 0's to make it 00n.  
        next_index = str("I"+(str(curr_index-1).zfill(3)))

        #loop back through list
        try:
            self.table.view.selection_set(next_index)
            self.table.view.focus(next_index)
            print(f"Moved focus: {next_index}")
        except:
            #go back to the bottom of the list and loop back
            max_row = str("I"+str(len(self.table.get_rows())).zfill(3))
            self.table.view.selection_set(max_row)
            self.table.view.focus(max_row)
            
            
if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()