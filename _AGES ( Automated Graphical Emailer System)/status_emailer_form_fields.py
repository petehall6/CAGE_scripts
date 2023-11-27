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

        
        self.pi = tbs.StringVar(value="")
        self.requester = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.objective = tbs.StringVar(value="")
        self.weeks = tbs.StringVar(value="")
        self.status_choice = tbs.StringVar(value="")
        
        self.create_labels()
        self.create_txtboxes()
        self.create_radiobtns()
        self.create_buttons()
        
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
        
        self.requester_lbl = tbs.Label(
            master = self.label_container, 
            text="Requested By: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.gene_lbl = tbs.Label(
            master = self.label_container, 
            text="Gene: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.cell_line_lbl = tbs.Label(
            master = self.label_container, 
            text="Cell Line: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.objective_lbl = tbs.Label(
            master = self.label_container, 
            text="Objective: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.weeks_lbl = tbs.Label(
            master = self.label_container, 
            text="Weeks: ", 
            font=(10), 
            bootstyle = SUCCESS,
        )
                
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.pi_lbl.grid(column=0,row=1, pady=20, sticky=W)
        self.requester_lbl.grid(column=0,row=2, pady=20, sticky=W)
        self.gene_lbl.grid(column=0,row=3, pady=20, sticky=W)
        self.cell_line_lbl.grid(column=0,row=4, pady=20, sticky=W)   
        self.objective_lbl.grid(column=0,row=5, pady=20, sticky=W)
        self.weeks_lbl.grid(column=0,row=6, pady=20, sticky=W)            
        
    def create_txtboxes(self):
        self.pi_box = tbs.Entry(
            master = self.label_container,
            bootstyle =PRIMARY,
            textvariable=self.pi,
        ) 
        self.requester_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
            textvariable=self.requester,
        )
        
        self.gene_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
            textvariable=self.gene,
        )
        
        self.cell_line_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
            textvariable=self.cell_line,
        )
        
        self.objective_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
            textvariable=self.objective,
        )
        
        self.weeks_box = tbs.Entry(
            master = self.label_container,
            bootstyle = PRIMARY,
            textvariable=self.weeks,
        )

        self.pi_box.grid(column=1, row=1, pady=10)
        self.requester_box.grid(column=1, row=2, pady=10)
        self.gene_box.grid(column=1, row=3, pady=10)
        self.cell_line_box.grid(column=1, row=4, pady=10)
        self.objective_box.grid(column=1, row=5, pady=10)
        self.weeks_box.grid(column=1, row=6, pady=10)

    def create_radiobtns(self):
        
        self.pool_radiobtn = tbs.Radiobutton(
            master = self.label_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Confirmed Pool",
            value = "Confirmed Pool",
        )
        
        self.screen_radiobtn = tbs.Radiobutton(
            master = self.label_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Initial Screen",
            value = "Initial Screen",
        )
        
        self.delay_radiobtn = tbs.Radiobutton(
            master = self.label_container,
            bootstyle = "info",
            variable = self.status_choice,
            text = "Delayed",
            value = "Delayed",
        )
        
        
        self.pool_radiobtn.grid(column=2,row=1,sticky=W)
        self.screen_radiobtn.grid(column=2,row=2,sticky=W)
        self.delay_radiobtn.grid(column=2,row=3,sticky=W)
        
    def create_buttons(self):

        self.gen_emails_btn = tbs.Button(
            master = self.label_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = PRIMARY,
            width = 25   
        )
        
        self.clear_btn = tbs.Button(
            master = self.label_container,
            text = "Clear Entries",
            command = self.clear_entries,
            bootstyle = DANGER,
            width = 25   
        )
        
        self.gen_emails_btn.grid(column=1, row=7, pady=40)
        self.clear_btn.grid(column=1, row=8, pady=40)

    def generate_emails(self):
        
        def _get_text_boxes(self):

            entry_input = [self.pi, self.requester, self.gene, self.cell_line, self.objective, self.weeks, self.status_choice]
            entries=[]
            
            
            for input in entry_input:
                if len(input.get()) != 0:
                    entries.append(input.get())
                    self.title_lbl.configure(
                        text = "Status Emailer",
                        font = ('Sans',25,'bold'),
                        bootstyle = WARNING,
                    )
                else:
                    self.title_lbl.configure(
                        text='Please Complete All Fields',
                        font = ('Sans',25,'bold'),
                        bootstyle=DANGER,
                    )
            
            pi,requester,gene,cell_line,objective,weeks,status = entries
            
            return pi, requester, gene, cell_line, objective, weeks, status 
        
        def _get_subject_line(gene, cell_line, objective):
            
            sub_line = f"{gene} {objective} {cell_line} status update"
                    
            return sub_line

        def _body_builder(status,greeting,cell_line,gene,objective,weeks):
            
            if status.upper() == "CONFIRMED POOL":

                body=f"""{greeting},
                <br><br>
                Great News! We have successfully confirmed the desired edit in the cell pool for your {cell_line} {gene} {objective} project.  
                We have already sorted for single cells into 96-well plates and will update you when we have screened 
                the plates and identified correctly edited clones. Each modification and cell line is a custom project, 
                and the time will differ widely for each project depending on several factors. Based on the details of your 
                specific project, we estimate that we will have identified clones in about {weeks} weeks. 
                If you have any questions or concerns, please don't hesitate to reach out.
                <br><br>
                """
                
            elif status.upper() == "INITIAL SCREEN":
                body=f"""{greeting},
                <br><br>
                Great News! We have successfully identified clones with the desired modification for your {cell_line} {gene} {objective} project. 
                If everything goes as planned, we expect to hand off these clones to you in {weeks} weeks. 
                Please let me know if you have any questions.
                <br><br>
                """
                
            elif status.upper() == "DELAYED": 
                body=f"""{greeting},
                <br><br>
                I wanted to provide you with an update on your {cell_line} {gene} {objective} project. Unfortunately, we were unable to identify any correctly 
                edited clones during the initial screen.  We are reviewing our data and reevaluating the editing strategy now.  
                We are still working hard to obtain the desired edited clone(s), but there is going to be a delay in the timeline 
                as we restart the process.  Please let me know if you have any questions.
                <br><br>
                """
                
            return body
                        
        signature = parse_signature()
        pi, requester, gene, cell_line, objective, weeks, status = _get_text_boxes(self)
        
        outlook = win32com.client.Dispatch("Outlook.Application")
        email = outlook.CreateItem(0)
        
        #removes duplicates and rephrases the greeting to a single person
        recip_list = [requester,pi]
        email_recip = list(set(recip_list))
        
        if len(email_recip) > 1:
            greeting = f"Hi {pi.split(' ')[0]} and {requester.split(' ')[0]}"
        else:
            greeting = f"Hi {pi.split(' ')[0]}"
        
        email_sub = _get_subject_line(gene,cell_line, objective)

        body = _body_builder(status,greeting,cell_line,gene,objective, weeks)

        email.To = ";".join(email_recip)
        email.CC = "Shondra Miller" #";".join(email_cc).replace(".","")

        email.bcc = "Shaina Porter"
        email.Subject = email_sub

        
        email.HTMLBody = body + signature
        #Display(False) loads all emails at once and gives focus back to ttk window
        email.Display(False)

    def clear_entries(self):
        entry_boxes = [self.pi_box, self.requester_box, self.gene_box, self.cell_line_box, self.objective_box, self.weeks_box]
        entry_vars = [self.pi, self.requester, self.gene, self.cell_line, self.objective, self.weeks, self.status_choice]
        
        for box in entry_boxes:
            box.delete(0, END)
        
        for var in entry_vars:
            var.set("")

        self.pi_box.focus_set()
    
if __name__ == '__main__':
    from _emailer_gui_RUN_THIS_SCRIPT import *
    app.mainloop()
    