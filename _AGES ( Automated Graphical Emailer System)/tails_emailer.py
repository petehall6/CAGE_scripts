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
import re
import os
import win32com.client
import glob
import datetime
from itertools import chain
import time

class Tails_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.header_container = tbs.Frame(self)
        self.header_container.pack(side=TOP,fill=X, expand=YES, pady=(15,10))
        
        self.table_container = tbs.Frame(self)
        self.table_container.pack(fill=BOTH,expand=YES, pady=5)
        
        self.button_container = tbs.Frame(self)
        self.button_container.pack(side=BOTTOM, fill=X, expand=YES, pady=(10,10))
        
        self.excel_name = tbs.StringVar(value=" ")
        self.srm_order = tbs.StringVar(value=" ")
        self.PI = tbs.StringVar(value=" ")
        self.requested_by = tbs.StringVar(value=" ")
        self.project_number = tbs.StringVar(value=" ")
        self.project_scope = tbs.StringVar(value=" ")
        self.cell_line = tbs.StringVar(value=" ")
        self.project_objective = tbs.StringVar(value=" ")
        self.gene = tbs.StringVar(value=" ")
        self.pi_department = tbs.StringVar(value=" ")
        self.success_choice = tbs.StringVar(value=" ")
        self.success_num = tbs.StringVar(value=" ")
        self.submitted_num = tbs.StringVar(value=" ")
        self.edit_choice = tbs.StringVar(value=" ")
        self.edit_size =tbs.StringVar(value=" ")
        self.injection_core = tbs.StringVar(value=" ")
        self.cage_nums = tbs.StringVar(value=" ")
        self.cage_dir = tbs.StringVar(value=" ")
        self.dir_selected = tbs.StringVar(value=" ")
        self.selected_programs = tbs.StringVar(value=" ")
    
        
        self.data = []
        self.cage_data = []
        self.selected_projects = []
        
        
        
        self.create_labels()
        self.create_buttons()
        self.create_comboxes()
        self.create_ngs_datepicker()

        self.table = self.create_table()
        self.cage_table = self.create_cage_table()

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
            text = 'Update',
            command =self.store_clicked,
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
        
        self.find_crispy_btn = tbs.Button(
            master=self.button_container,
            text='Find CRISPY Programs',
            command=self.find_crispy_files,
            bootstyle = SECONDARY,
            width = 20
        )
        
        self.confirm_crispy_btn = tbs.Button(
            master = self.button_container,
            text="Confirm CRISPY Programs",
            command = self.confirm_crispy_click,
            bootstyle = SUCCESS,
            width = 25
        )
        
        self.srm_load_btn.grid(column=0,row=1, pady=10)
        self.store_btn.grid(column=1, row=12,pady=10,sticky=W)
        self.prev_btn.grid(column=0, row=13,pady=10,padx=5)
        self.next_btn.grid(column=2, row=13,pady=10,padx=5)
        self.gen_emails_btn.grid(column=0, row=14, pady=60)
        self.clear_btn.grid(column=2, row=14,pady=60, sticky=N+S+E+W, padx = 10)
        self.find_crispy_btn.grid(column=6, row=6, pady=10,padx=10, sticky=W)
        self.confirm_crispy_btn.grid(column=6, row=6, pady=10, padx=10, sticky=E)
    
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
        
        self.cage_num_lbl = tbs.Label(
            master = self.button_container,
            text = 'Enter CAGE Numbers:',
            font=(10),
            bootstyle = SUCCESS,
        )
        
        self.ngs_date_lbl = tbs.Label(
            master = self.button_container,
            text = 'Select NGS Date:',
            font=(10),
            bootstyle = SUCCESS,
        )
        
        self.ngs_date_error_lbl = tbs.Label(
            master = self.button_container,
            text="",
            font = (8),
            bootstyle = DANGER
            
        )
        
        self.crispy_status_lbl = tbs.Label(
            master = self.button_container,
            text=" ",
            font=(10),
            bootstyle = INFO
        )
        
        self.success_num_lbl = tbs.Label(
            master = self.button_container,
            text="Num. Succeeded",
            font=(10),
            bootstyle = SUCCESS
        )
        
        self.submitted_num_lbl = tbs.Label(
            master = self.button_container,
            text="Num. Submitted:",
            font=(8),
            bootstyle = SUCCESS
        )
        
        self.edit_size_lbl = tbs.Label(
            master = self.button_container,
            text="Size of Edit:",
            font=(8),
            bootstyle = SUCCESS
        )
        
        self.injection_lbl = tbs.Label(
            master = self.button_container,
            text="TCU vs NEL:",
            font=(10),
            bootstyle = SUCCESS
        ) 
        
        self.pi_department_lbl = tbs.Label(
            master = self.button_container,
            text = "Department: ",
            font = (10),
            bootstyle = SUCCESS
        )
        
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10, sticky=W)
        
        self.proj_lbl.grid(column = 0,row=3, sticky=E)
        self.active_proj_lbl.grid(column=1,row=3, columnspan=3, sticky=W+E+N+S)
        
        self.success_lbl.grid(column=0,row=5, pady=10)
        self.edit_lbl.grid(column=0, row=7, pady=10)
        self.cage_num_lbl.grid(column=0,row=8, pady=10)
        self.ngs_date_lbl.grid(column=0, row=9, pady=10)
        self.ngs_date_error_lbl.grid(column=2, row=9)
        self.crispy_status_lbl.grid(column=6, row=7,sticky=W+E)
        self.injection_lbl.grid(column=0, row=10, pady=10)
        self.pi_department_lbl.grid(column=0, row=11, pady=10)
        self.success_num_lbl.grid(column=2,row=4,sticky=S)
        self.edit_size_lbl.grid(column=2,row=6,sticky=S)
        self.submitted_num_lbl.grid(column=3,row=4,sticky=S)
        
    def create_comboxes(self):
        
        outcomes = [' ','Yes','No']
        edits = [' ','KO','KI','CKO','Del','ssODN','PM','Data']
        injection = [' ', 'TCU', 'NEL']
        department = ['CBT', 'CMB', 'Comp Bio', 'DNB', 'Pharm Sciences', 'Struct Bio', 'Tumor Bio', 'BMT', 'Infect Dis' ]
        
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
        
        self.cage_box = tbs.Entry(
            master = self.button_container,
            bootstyle = "info",
            width=25    
        )
            
        self.success_num_box = tbs.Entry(
            master = self.button_container,
            bootstyle = "info",
            width=25    
        )

        self.submitted_num_box = tbs.Entry(
            master = self.button_container,
            bootstyle = "info",
            width=25    
        )
        
        
        self.edit_size_box = tbs.Entry(
            master = self.button_container,
            bootstyle = "info",
            width=25    
        )

        self.injection_box = tbs.Combobox(
            master = self.button_container,
            bootstyle = "info",
            value = injection,
        )
        
        self.pi_department_box = tbs.Combobox(
            master = self.button_container,
            bootstyle = "info",
            value = department,
        )
        
        self.success_box.grid(column=1,row=5, sticky=W+E)
        self.edits_box.grid(column=1, row=7, sticky=W+E)
        self.cage_box.grid(column=1, row=8, sticky=W+E)
        self.success_num_box.grid(column=2, row=5)
        self.submitted_num_box.grid(column=3, row=5)
        self.edit_size_box.grid(column=2, row=7)
        self.injection_box.grid(column=1, row=10)
        self.pi_department_box.grid(column=1, row=11)
            
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
            {"text": 'Department'}, 
            {"text":'Success'},
            {"text":'Num Success'},
            {"text":'Num Submitted'},
            {"text":'Edit'},
            {"text":'Edit Size'},
            {"text":'TCU/NEL'},
            {"text":'Selected Programs'}
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

    def create_cage_table(self):
        
        columns = [
            {"text":'Programs'},
            {"text":'Selected'}
        ]
        
        self.cage_table = Tableview(
            master = self.button_container,
            coldata=columns,
            rowdata=self.cage_data,
            paginated=False,
            searchable=False,
            bootstyle=PRIMARY,
            stripecolor=LIGHT,
            autoalign=False,
        )
    
        #self.table.view.selection_set(0)
        self.cage_table.view.bind()
        self.cage_table.view.bind("<<TreeviewOpen>>", self.cage_table_clicked)
        self.cage_table.grid(column=6,row=3, rowspan=3)
        
        return self.cage_table

    def create_ngs_datepicker(self):
        
        self.ngs_date_picker = tbs.DateEntry(
            master = self.button_container,
            dateformat="%m%d%y",
        )

        self.ngs_date_picker.grid(column=1, row=9)
  
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
            #* its a tuple dummy                   
            self.data.append((self.srm_order,
                            self.PI,
                            self.requested_by,
                            self.entered_by,
                            self.project_number,
                            self.gene,
                            self.sample_num,
                            self.sample_format,
                            self.sample_type,
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
        
        #how to check if project has been editied in gui
        if len(selected_proj_info) != 9:
            print("setting combo value")
            self.success_box.set(selected_proj_info[9])
            self.edits_box.set(selected_proj_info[10])
        else:
            self.success_box.set(" ")
            self.edits_box.set(" ")

    def cage_table_clicked(self,event):
        
        program_select = self.cage_table.view.item(self.cage_table.view.focus(), "values"[0])
        
        if program_select[1] == ' ':
            pick_update = "Yes"
        elif program_select[1] == "Yes":
            pick_update = "No"
        elif program_select[1] == "No":
            pick_update = "Yes"
        
        
        #get row tuple where focus is (highlighted row)
        selected_proj_info = self.cage_table.view.item(self.cage_table.view.focus(),"values")
        
        proj_appended=[]
        
        
        for info in selected_proj_info:
                proj_appended.append(info)
        
        proj_appended[1] = pick_update
        
        #convert back to tuple to plug back into item()
        updated_proj_info = tuple(proj_appended)
        #text is normally blank, have to specifiy value update
        self.cage_table.view.item(self.cage_table.view.focus(),text="",values=updated_proj_info)

    def generate_emails(self):
        
        signature = parse_signature()
        
        def _email_writer(project_details,initial_choice):
            
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
                    #convert back into list and pass to writer
                    #kept in list form since single mode was originally written and multi project was an added on feature
                    proj_details = pi_specifc_df.values.tolist()[0]
                    
                    _email_writer(proj_details)
                    
        return

    def clear_controls(self):
        
        table = self.table
        
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        self.cage_table.delete_rows()
        
        self.success_box.set("")
        self.edits_box.set(" ")
        self.cage_box.delete(0, 'end')
        
        
        self.cage_table.unload_table_data()
        self.table.unload_table_data()
        self.data = []
        self.cage_data = []

    def store_clicked(self):
        
        #get combobox values
        self.pi_department = self.pi_department_box.get()
        self.success_choice = self.success_box.get()
        self.edit_choice = self.edits_box.get()
        self.edit_size = self.edit_size_box.get()
        self.success_num = self.success_num_box.get()
        self.submitted_num = self.submitted_num_box.get()
        self.injection_core = self.injection_box.get()
        #returns string. will need to split at comma and append CAGE where necessary
        self.cage_nums = self.cage_box.get().strip()
        
        cage_num_list = self.cage_nums.split(",")
        
        #check if any letters are in the item and replace with corrected format
        for project in cage_num_list:
            if re.search('[a-zA-Z]', project) == None and len(cage_num_list)>=1:
                cage_num_list = list(map(lambda x: x.replace(project,"CAGE"+project), cage_num_list))
        
        cage_nums_str = ",".join(cage_num_list)
        
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
        
        proj_appended[4] = cage_nums_str
        proj_appended.append(self.success_choice)
        proj_appended.append(self.success_num)
        proj_appended.append(self.submitted_num)
        proj_appended.append(self.edit_choice)
        proj_appended.append(self.edit_size)
        proj_appended.append(self.injection_core)
        
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
            
        #check tableview values and populate crispy table if there
        
        focused_project_files = self.table.view.item(self.table.view.focus(),"values")
        
        print(f"next button: {focused_project_files}")
        
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
    
    def find_crispy_files(self):
        
        self.cage_table.unload_table_data()
        self.cage_data=[]
        cage_dirs = []
        
        self.ngs_date_error_lbl.configure(text="")
        
        NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"
        
        
        ngs_run_date = self.ngs_date_picker.entry.get()
        
        #go to NGS date
        try:       
            os.chdir(os.path.join(NGS_DIR,ngs_run_date,"joined"))
            self.ngs_date_error_lbl.configure(text="NGS run found.", bootstyle=SUCCESS)        
            print(os.getcwd())
        except:
            self.ngs_date_error_lbl.configure(text="NGS run not found. Pick another date.")
        
        self.crispy_status_lbl.configure(text="Searching...")
        #get focus of tableview and get values tuples
        selected_proj_info = self.table.view.item(self.table.view.focus(),"values"[0])
        #parse CAGE numbers
        cage_nums = selected_proj_info[4].split(",")
        print(cage_nums)
        
        for num in cage_nums:
            self.crispy_status_lbl.configure(text=f"Searching {num}")
            #find all matches with CAGe# and only return directories
            cage_dirs.append(list([name for name in glob.glob(f"{num}*") if os.path.isdir(name)]))

        print(f"cage dirs: {cage_dirs}")

        #Add cage_folders to the cage_table
        
        #append empty index to list for selected spot

        #check if sub_list is empty ie didnt return hits
        if len(cage_dirs[0]) > 0:
            for proj in cage_dirs:
                #add in the empty string for the selected column 
                dir_list = list(proj)
                for item in dir_list:
                    cage_tup = tuple([item," "])
                #TODO
                    print(f"dir_list{dir_list}")
                    print(f"cage_tup: {cage_tup}")
                #and convert back to tuple since treeview items must be tuples
                    self.cage_data.append(cage_tup)
            self.cage_table.destroy()
            self.cage_table.load_table_data()
            self.cage_table = self.create_cage_table()
            self.crispy_status_lbl.configure(text="Project Folders Found")
            
        else:
            self.cage_table.unload_table_data()
            self.cage_table.destroy()
            self.create_cage_table()
            self.crispy_status_lbl.configure(text="No matching projects found.")
        
    def confirm_crispy_click(self):
        proj_appended=[]
        selected_cage_proj_info = self.table.view.item(self.table.view.focus(),"values")
        for info in selected_cage_proj_info:
            proj_appended.append(info)
        
        print(f"Len of proj_appended: {len(proj_appended)}")
        #get info from cage_table if rows are marked as Yes
        cage_programs_info = self.cage_table.view.get_children()
        found_programs=[]
        selected_cage_program=[]
        
        print(f"proj_appended {proj_appended}")
        
        #check if programs have been added:
        if len(proj_appended) == 16:
            for item in cage_programs_info:
                if self.cage_table.view.item(item)['values'][1] == 'Yes':
                    selected_cage_program.append(self.cage_table.view.item(item)['values'][0])
                    proj_appended.append(selected_cage_program)
                print(f"You chose me!!: {selected_cage_program}")

            #convert back to tuple to plug back into item()
            updated_proj_info = tuple(proj_appended)
            
        else:
            selected_cage_program=[]
            proj_appended[-1] = []
            for item in cage_programs_info:
                if self.cage_table.view.item(item)['values'][1] == 'Yes':
                    selected_cage_program.append(self.cage_table.view.item(item)['values'][0])
                    proj_appended[-1] = selected_cage_program
                print(f"You chose me!!: {selected_cage_program}")
            
        updated_proj_info = tuple(proj_appended)
            
        #text is normally blank, have to specifiy value update
        self.table.view.item(self.table.view.focus(),text="",values=updated_proj_info)

        
        print(f"updated: {updated_proj_info}")
    
    
    def generate_emails(self):
        
        signature = parse_signature()
        
        def _email_writer(project_details):
            
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
                    Your live cultures will be on the bottom shelf of the "Pick-up" incubator, which is labeled accordingly.  Please bring dry ice for the pickup.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
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
                    Your live cultures will be on the bottom shelf of the "Pick-up" incubator, which is labeled accordingly.  Please bring dry ice for the pickup.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>
                    """
                    
                return body
            
            def _update_mice_excel():
                
                
                return
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

                    #print(f"single project: {pi_specifc_df.iloc[0][1]}")
                    #convert back into list and pass to writer
                    #kept in list form since single mode was originally written and multi project was an added on feature
            proj_details = pi_specifc_df.values.tolist()[0]
            
            _email_writer(proj_details)
                    
        return
   
        
if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()