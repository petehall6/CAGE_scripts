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


class DropOff_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.project_scope = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.project_objective = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.species = tbs.StringVar(value="")
        self.line_lead = tbs.StringVar(value="")
        self.stem_cell = tbs.StringVar(value="")
        
        
        self.data = []
        
        self.create_srm_load_btn()
        self.create_labels()

        
        self.table = self.create_table()
        self.create_gen_emails_btn()
        self.create_clear_btn()

        
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
            {"text": "Species"},
            {"text":'Project Scope'},
            {"text":'Cell Line of Choice'},
            {"text": 'Project Objective'},
            {"text":'Target Gene Name'},
            {"text":'Cell Line Lead'},
            {"text": 'Stem Cells?'}
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
                self.species = srm[8]
                self.line_lead = srm[9].split("Lead-")[1]
                self.stem_cell = srm[10]
                
            self.data.append((self.srm_order,
                              self.PI,
                              self.requested_by,
                              self.project_number,
                              self.species,
                              self.project_scope,
                              self.cell_line,
                              self.project_objective,
                              self.gene,
                              self.line_lead,
                              self.stem_cell))
            
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        def _get_subject_line(scope,species,gene,cell_line, objective, stem_cell):
            
            if scope.upper() == "CELL LINE CREATION" and species.upper() == "HUMAN" and stem_cell.upper() == "NO":
                sub_line = f"{species} {cell_line} cell line intake"
            elif scope.upper() == "CELL LINE CREATION" and species.upper() == "MOUSE" and stem_cell.upper() == "NO":
                sub_line = f"{species} {cell_line} cell line intake"
            elif stem_cell.upper() == YES:
                sub_line = f"{cell_line} cell line intake"
                
            return sub_line

        def _body_builder(requester, pi, scope, cell_line, objective, line_lead, stem_cell):

            pi = pi.split(", ")[1]
            requester = requester.split(", ")[1]            
            
            if scope.upper() == "CELL LINE CREATION" and species.upper()  == "HUMAN" and stem_cell.upper() == "NO":
                body=f"""Hi {requester},
                <br><br>
                We are ready to intake the {cell_line} cells for your {gene} {objective} projects.
                <br><br>
                We have a contactless drop off system in place.  
                Please arrange for someone to drop off the items below in the new ARC building, 4th floor, M4170.  
                To find the CAGE, take the elevators to the 4th floor, and turn right at the first two hallways. 
                We are at the end of the hallway. 
                The live cells can go in our quarantine incubator, which can be found in the right side of M4170, on the floor under the shelves before the hoods. 
                It is labeled as quarantine incubator. The door is always unlocked.  The media can go in the same room, in the labeled fridge to the right of the quarantine incubator.  
                If you need help, feel free to ask anyone in the CAGE.
                <br><br> 
                Once drop off is complete, please email {line_lead.split(" ")[0]} to let them know.
                <br><br>
                1.	T75 flask of live cells<br>
                2.	500 ml of complete media<br>
                3.	An electronic copy of the media recipe and any special culturing conditions<br>
                4.	A recent (within the last 3 months) STR profile from the Hartwell Center<br>
                <br><br>
                Thanks,
                <br><br>
                Shaina               
                """
                
            elif scope.upper() == "CELL LINE CREATION" and species.upper() == "MOUSE" and stem_cell.upper() == "NO":
                body=f"""Hi {requester},
                <br><br>
                We are ready to intake the {cell_line} cells for your {gene} {objective} projects.
                <br><br>
                We have a contactless drop off system in place.  
                Please arrange for someone to drop off the items below in the new ARC building, 4th floor, M4170.  
                To find the CAGE, take the elevators to the 4th floor, and turn right at the first two hallways. 
                We are at the end of the hallway. 
                The live cells can go in our quarantine incubator, which can be found in the right side of M4170, on the floor under the shelves before the hoods. 
                It is labeled as quarantine incubator. The door is always unlocked.  The media can go in the same room, in the labeled fridge to the right of the quarantine incubator.  
                If you need help, feel free to ask anyone in the CAGE.
                <br><br>
                Once drop off is complete, please email {line_lead.split(" ")[0]} to let them know.
                <br><br>
                1.	T75 flask of live cells<br>
                2.	500 ml of complete media<br>
                3.	An electronic copy of the media recipe and any special culturing conditions<br>
                <br><br>
                Thanks,
                <br><br>
                Shaina
                """
                
            elif stem_cell.upper() == "YES": 
                body=f"""We are ready to intake the {cell_line} cells for your {gene} {objective} project.
                <br><br> 
                We have a contactless drop off system in place.  
                Please arrange for someone to drop off the items below in the new ARC building, 4th floor, M4170.  
                To find the CAGE, take the elevators to the 4th floor, and turn right at the first two hallways.
                We are at the end of the hallway. 
                The live cells can go in our quarantine incubator, which can be found in the right side of M4170, on the floor under the shelves before the hoods.
                It is labeled as quarantine incubator. The door is always unlocked.  If you need help, feel free to ask anyone in the CAGE. 
                <br><br>
                Once drop off is complete, please email {line_lead.split(" ")[0]} to let them know.
                <br><br>
                1.	6 wells of a 6 well plate of live cells, 40% confluent<br>
                2.	An electronic copy of any special culturing conditions<br>
                3.	A recent (within the last 3 months) STR profile from the Hartwell Center<br>
                4.	A recent (within two months of the freeze date) karyotype<br>
                <br><br>
                Thanks,
                <br><br>
                Shaina
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
            

            srm_order_num, pi, requester, project_num,species, scope, cell_line, objective, gene, line_lead, stem_cell = entry

            #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            email_recip = str(requester)
            email_cc = str(pi)
            
            email_sub = _get_subject_line(scope,species,gene,cell_line, objective, stem_cell)

            body = _body_builder(requester,pi,scope,cell_line,objective, line_lead, stem_cell)

            email.To = email_recip
            email.CC = email_cc

            #email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)


if __name__ == "__main__":
    import emailer_gui
    emailer_gui.app.mainloop()