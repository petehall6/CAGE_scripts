import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from tkinter import ttk

from emailer_functions import (open_file,
                               df_from_template,
                               parse_signature
                            )
import pandas as pd
import os
import win32com.client
import glob



class Batch_Tab(tbs.Frame):
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
        self.email_template = tbs.StringVar(value="")

        style = ttk.Style()
        style.configure("Custom.Treeview",rowheight=30)

        
        
        self.data = []

        self.create_labels()
        self.create_srm_load_btn()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        self.table = self.create_table()

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Batch Emailer",
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
            {"text":"Template"},          
            {"text":'Gene'},
            {"text":'gRNA 1'},
            {"text":'gRNA 2'},
            {"text":'gRNA 3'},
            {"text":'gRNA 4'},
            {"text":'PD name'},
            {"text":'ssODN#1 name'},
            {"text":'ssODN#2 name'}
        ]

        self.table = Tableview(
            master = self,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=False,
            bootstyle=PRIMARY,
            stripecolor=DARK,)
        
        
        self.table.view.configure(style="Custom.Treeview")
        self.table.view.bind("<ButtonRelease-1>", self.on_select)
        self.table.pack(side=BOTTOM,fill=BOTH, expand=YES, padx=10, pady=10)

        return self.table

    def clear_controls(self):
        
        table = self.table
        #clear label and loaded template
        self.excel_name = tbs.StringVar(value="")
        self.excel_lbl.config(text="")
        
        table.delete_rows()
        widgets = self.table.winfo_children()
        print(widgets)
        for widget in widgets:
            print(widget.winfo_class())
            if widget.winfo_class() == 'TCombobox':
                widget.destroy()
        
        
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
        
        
        self.email_template = "test"
        
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
                #self.email_template will load from dropdown list
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
                            self.email_template,
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

    def on_select(self, event):
        print("clicked")
        item_id = self.table.view.selection()[0]
        column = self.table.view.identify_column(event.x)
        print(f"column: {column}.  item_id: {item_id}")
        print(event)
        print(event.x)
        values =[
            "gRNA with validation", 
            "gRNA and donor - gRNA validation",
            "gRNA and donor - donor validation", 
            "Animal model - gRNA GEMM",
            "Animal model - gRNA NEL",
            "Animal model - AAV GEMM",
            "Animal model - cKO GEMM",
            "Animal model - donor NEL",
        ]
        if column == "#7":  # Check if the click is in the first column
            x, y, width, height = self.table.view.bbox(item_id, column)
            combo = tbs.Combobox(self.table, values=values)
            combo.place(x=x, y=y, width=width, height=height)
            combo.current(0)  # Set the initial value
            combo.bind("<<ComboboxSelected>>", lambda e: self.update_cell(item_id,column, combo.get()))
            combo.focus_set()
            
            print("Clicked in the column")

    def update_cell(self,item_id,column,value):
        current_values = list(self.table.view.item(item_id, "values"))
        col_index = int(column.replace("#","")) - 1
        current_values[col_index] = value
        
        print(f"current_values: {current_values}")
        
        self.table.view.item(item_id, values=current_values)
        widgets = self.table.winfo_children()
        
        self.data = self.table.get_rows(visible=True)
        
        row_index = self.table.view.index(item_id)
        self.data[row_index] = current_values
        
        
        print(f"Updated self.data: {self.data}")

        for widget in widgets:
            if widget.winfo_class() == 'TCombobox':
                widget.destroy()
        pass        

    def generate_emails(self):
        def _get_subject_line(scope):
            
                if scope.lower() == "animal model creation":
                    sub_line = f"test__animal_model_creation__test"
                    
                elif scope.lower() == "grna and donor with validation":
                    sub_line = f"test__grna_and_donor_w_val__test"

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
        
        def _body_builder(greeting, objective, template):
            if template.lower() == "grna with validation":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  We have identified active CRISPR nucleases for your {objective} project.  Please see the updated slide deck.  
                <br><br>
                I am packaging your reagents, and they will be ready for pick up anytime after tomorrow.
                <br><br>
                When you walk into M4160, there are two mini-freezers under the bench on the left. In the door of the freezer on the right, there is a clear accordion folder, and in that folder your guides will be in a baggie under the first letter of your PI's last name.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                Best,
                <br><br>
                <br><br>
                </font>                
                """
                
            elif template.lower() == "grna and donor - grna validation" or template.lower() == "animal model - gRNA GEMM" or template.lower() =="animal model - gRNA NEL":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  We have identified active CRISPR nucleases for your {objective} project.  Please see the updated slide deck.  We will now move on to donor design and validation.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                Thanks!
                <br><br>
                <br><br>
                </font>
                """
                
            elif template.lower() == "grna and donor - donor validation": 
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  We have completed donor validation for your {objective} project.  Please see the updated slide deck.  
                <br><br>
                I am packaging your reagents, and they will be ready for pick up anytime after tomorrow.
                <br><br>
                When you walk into M4160, there are two mini-freezers under the bench on the left. In the door of the freezer on the right, there is a clear accordion folder, and in that folder your guides will be in a baggie under the first letter of your PI's last name.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                Thanks!
                <br><br>
                <br><br>
                </font>
                """
            
            elif template.lower() == "animal model - aav gemm":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  We have completed donor validation for your {objective} project.  Please see the updated slide deck. 
                <br><br>
                To increase efficiency, we will now schedule zygote manipulations with GEMM on your behalf using the cost center associated with the initial CAGE request.  If you do not wish to move forward with creating this mouse model, please let us know in the next 72 hours.  Otherwise, we will move forward with mouse model creation.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                <br><br>
                Thanks!
                </font>
                """
            
            elif template.lower() == "animal model - cko gemm":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  Your donor AAV particles are ready for your {objective} project, and we are ready to move forward with animal model creation.
                <br><br>
                To increase efficiency, we will now schedule zygote manipulations with GEMM on your behalf using the cost center associated with the initial CAGE request.  If you do not wish to move forward with creating this mouse model, please let us know in the next 72 hours.  Otherwise, we will move forward with mouse model creation.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                Thanks!
                <br><br>
                <br><br>
                """
            elif template.lower() == "animal model - donor nel":
                body=f"""
                <font face="Calibri, Calibri, monospace">
                {greeting},
                <br><br>
                Great news!  We have completed donor validation for your XXX project.  Please see the updated slide deck. 
                <br><br>
                When you are ready, please reach out to Valerie Stewart to schedule your zygote manipulations.
                <br><br>
                Please let us know if you have any questions.
                <br><br>
                Thanks!
                <br><br>
                <br><br>
                """
            else:
                body =  "Double check project scopes and objectives are correct in the SRM Excel Export file."
            return body

        signature = parse_signature()
        
        intact_rows = self.table.get_rows(visible=True)
        table_entries=[]
        for row in intact_rows:
            table_entries.append(row)
        
            """
            template options:
        
            "gRNA with validation", 
            "gRNA and donor - gRNA validation",
            "gRNA and donor - donor validation", 
            "Animal model - gRNA GEMM",
            "Animal model - gRNA NEL",
            "Animal model - AAV GEMM",
            "Animal model - cKO GEMM",
            "Animal model - donor NEL",

            """
        
        for entry in table_entries:
            project_num, pi, requester, scope, species, objective, template, gene, grna_1, grna_2, grna_3, grna_4, pd_name, ssodn_1, ssodn_2 = entry
            
            objective = str(objective).lower()
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
            
            email_cc = []
                            
            email_sub = _get_subject_line(scope)

            #email = _get_attachment(email,project_num)

            body = _body_builder(greeting, objective, template)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc).replace(".","")

            email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            
            email.HTMLBody = body + signature
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)
                            
            return


if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()