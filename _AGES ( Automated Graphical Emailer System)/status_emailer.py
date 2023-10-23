import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (clicked,
                               parse_signature
                        )
import os
import win32com.client

class Status_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.label_container = tbs.Frame(self)
        self.label_container.pack(fill=BOTH, expand=YES)

        
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.status = tbs.StringVar(value="")
        
        self.status_choice = ['Confirmed Pool', 'Inital Screen', 'Delayed']
        
        self.create_labels()
        self.create_txtboxes()
        self.create_radiobtns()
        self.create_gen_emails_btn()
        
    def create_labels(self):
        self.title_lbl = tbs.Label(
            master = self.label_container,
            text = "Status Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
        
        self.pi_lbl = tbs.Label(
            master = self.label_container, 
            text="PI: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.requestedby_lbl = tbs.Label(
            master = self.label_container, 
            text="Requested By: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.pi_lbl.grid(column=0,row=1, pady=10)
        self.requestedby_lbl.grid(column=0,row=2, pady=10)      
        
    def create_txtboxes(self):
        self.pi_box = tbs.Entry(
            master = self.label_container,
            bootstyle =PRIMARY,
        ) 
        
        self.requestedby_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
        )

        self.pi_box.grid(column=1,row=1)
        self.requestedby_box.grid(column=1, row=2)

    def create_radiobtns(self):
        row_num=1
        for status in self.status_choice:
            tbs.Radiobutton(
                master = self.label_container,
                bootstyle = "info",
                variable = self.status_choice,
                text = status,
                value = status,
                command = ""
            ).grid(column=3,row=row_num, pady=5)
            row_num+=1
        
        
    def create_gen_emails_btn(self):

        self.gen_emails_btn = tbs.Button(
            master = self.label_container,
            text = "Create Emails",
            #command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.gen_emails_btn.grid(column=1, row=4, pady=30)


    def generate_emails(self):
        
        def _get_subject_line(scope, gene, cell_line, objective):
            
            if scope.upper() == "EDITED CELL POOL":
                sub_line = f"{gene} {cell_line} Edited Cell Pool Complete"
            elif scope.upper() == "CELL LINE CREATION":
                sub_line = f"{cell_line} {gene} {objective} Cell Line Complete"
            elif scope.upper() == "CELL FITNESS/DEPENDENCY ASSAY":
                sub_line = f"CelFi Assay for {gene} in {cell_line} Cells Complete"
                    
            return sub_line
        
        
        def _body_builder(requester, pi, scope, cell_line, objective, line_lead):
            
            pi = pi.split(", ")[1]
            requester = requester.split(", ")[1]
            
            
            if scope.lower() == "edited cell pool":

                body=f"""Hi {pi} and {requester},
                <br><br>
                Great news! Your {gene} {cell_line} edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.
                <br><br>
                The last slide is the most informative.  We were able to get over <font color=red>XX%</font> total editing in the pool with <font color=red>~XX%</font> out of frame indels.
                <br><br>
                We have a contactless pickup system in place right now.  Please coordinate with {line_lead.split(" ")[0]} to let them know a good time window for you to pick up these cells. 
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
                <br><br>
                """
                
            else: 
                body=f"""Hi {pi} and {requester},
                <br><br>
                Great news! Your {cell_line} {gene} {objective} project is complete and ready for pick up.  Please see the attached slide deck for details.
                <br><br>
                We currently have a contactless pickup system in place.  Please arrange a time window with {line_lead.split(" ")[0]} in which someone can pick up the cells.  
                At the agreed upon time, {line_lead.split(" ")[0]} will place your frozen vials of cells into a dry ice bucket in M4170.  
                The dry ice bucket will be on the counter in front of you as you walk in.  
                Your live cultures will be in the first incubator to the right (top incubator, bottom shelf) and labeled accordingly. Please also bring dry ice for the pickup.
                <br><br>
                As always, please let me know if you have any questions.
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
            email_recip = [requester]
            email_cc = [pi,line_lead]
            
            email_sub = _get_subject_line(scope,gene,cell_line, objective)

            body = _body_builder(requester,pi,scope,cell_line,objective, line_lead)

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc).replace(".","")

            email.bcc = "Shaina Porter"
            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)

if __name__ == '__main__':
    from _emailer_gui_RUN_THIS_SCRIPT import *
    app.mainloop()
    