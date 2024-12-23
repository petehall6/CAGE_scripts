import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (open_file,
                               df_from_design_template,
                               parse_signature
                            )
import pandas as pd
import os
import win32com.client
import glob


class Design_Tab(tbs.Frame):
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
        self.initial_choice = tbs.StringVar(value="")
        
        self.data = []

        self.create_labels()
        self.create_buttons()
        self.create_radiobtns()

        self.table = self.create_table()

    def create_buttons(self):
        self.srm_load_btn = tbs.Button(
            master = self.button_container,
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
        
        self.srm_load_btn.grid(column=0,row=1, pady=10)
        self.gen_emails_btn.grid(column=0, row=2, pady=10)
        self.clear_btn.grid(column=0, row=4,pady=60)
        
    def create_radiobtns(self):
        self.jk_radiobtn = tbs.Radiobutton(
        master = self.button_container,
        bootstyle = "info",
        variable = self.initial_choice,
        text = "JK",
        value = "JK",
    )
    
        self.bh_radiobtn = tbs.Radiobutton(
        master = self.button_container,
        bootstyle = "info",
        variable = self.initial_choice,
        text = "BH",
        value = "BH",
    )
    
        self.jk_radiobtn.grid(column=3,row=1,sticky=W)
        self.bh_radiobtn.grid(column=3,row=2,sticky=W)
        
    def create_labels(self):
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Designs Emailer",
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
        self.initial_choice.set("")

    def load_srm(self):
        self.table.unload_table_data()
        self.data=[]
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        srm_list = df_from_design_template(template)
        
        #append table to bottom of frame
        #loop through index
        #then loop through element
        
        try:
            for srm in srm_list:
                for entry in srm:
                    self.srm_order = srm[0]
                    self.project_number = srm[1].strip().split(" and")[0]
                    self.PI = srm[2]
                    self.requested_by = srm[3]
                    self.project_scope = srm[4]
                    self.gene = srm[5]
                    self.grna_num = srm[6]
                    self.grna1 = srm[7]
                    self.grna2 = srm[8]
                    self.grna3 = srm[9]
                    self.grna_final = srm[10]                    
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
            print(f"Project: {self.project_number}")
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
                Attached are the designs, off-target analysis, and our recommendations for which gRNAs to move forward 
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
                email_cc.append("Baranda Hansen")
            else:
                email_cc.append("Jon Klein")
            
            
            if len(email_recip) > 1:
                greeting = f"Hi {pi.split(',')[1]} and {requester.split(',')[1]}"
            else:
                greeting = f"Hi {pi.split(',')[1]}"
                            
            email_sub = _get_subject_line(gene,srm_order_num)

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
            
            def _get_subject_line(gene,srm_order_num):
                
                
                srm_order_str = ""
                
                for order in srm_order_num:
                    srm_order_str = srm_order_str + str(order) + ', '
            
                sub_line = f"gRNA Designs for {(', '.join(gene))}, SRM orders: {srm_order_str}"
                
                return sub_line
            
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
                
                body = f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Attached are the designs, off-target analysis, 
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
                </font>
                """

                return body
            
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = requester + pi
            
            email_recip = list(set(recip_list))
            
            if len(email_recip) > 1:
                first_names=[]
                
                for inv in pi:
                    first_names.append(inv.split(',')[1])
                
                for req in requester:
                    first_names.append(req.split(',')[1])
                first_name_list = list(set(first_names))
    
                #insert 'and' in front of last element
                            
                first_name_list[-1] = ' and '+ first_name_list[-1] 
                
                first_name_str =''
                
                for name in first_name_list:
                    first_name_str = first_name_str + str(name)
                    
                
                
                greeting = f"Hi {first_name_str}"


            else:
                greeting = f"Hi {pi[0].split(',')[1]}"
            
            #copies either Jon or Barnanda depending on who sent the email
            email_cc = ["Shondra Miller"]
            if initial_choice == "JK":
                email_cc.append("Baranda Hansen")
            else:
                email_cc.append("Jon Klein")
                            
            email_sub = _get_subject_line(gene,srm_order_num)

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
        

if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()