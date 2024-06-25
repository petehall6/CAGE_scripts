import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from ttkbootstrap.dialogs import Messagebox
from emailer_functions import (open_file,
                               df_from_ngs_template,
                               parse_signature
                            )
import os
import win32com.client
import glob


class NGS_Tab(tbs.Frame):
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
        self.gene = tbs.StringVar(value="")
        self.user_comments = tbs.StringVar(value="")
        self.ngs_date = tbs.StringVar(value="")

        
        
        self.data = []

        self.create_labels()
        
        self.create_srm_load_btn()
        self.create_ngs_date_field()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        
        self.table = self.create_table()

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "NGS Emailer",
            font = ('Sans',25,'bold'),
            bootstyle = WARNING,
        )
        
        
        self.excel_lbl = tbs.Label(
            master = self.button_container, 
            text="SRM template", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.ngs_lbl = tbs.Label(
            master = self.button_container, 
            text="NGS Date:", 
            font=(10), 
            bootstyle = SUCCESS,
        )
        
        self.title_lbl.grid(column=1,row=0, columnspan=3, padx=20, sticky=W+E+N+S)
        self.excel_lbl.grid(column=1,row=1, pady=10)
        self.ngs_lbl.grid(column=0, row=2)
    
    def create_ngs_date_field(self):
        
        self.ngs_date_box = tbs.Entry(
        master = self.button_container,
        bootstyle = "info",
        width=22
        )
            
        self.ngs_date_box.grid(column=1, row=2, sticky='w')
        
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
        
        self.gen_emails_btn.grid(column=0, row=3, pady=10)

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
            {"text":'Requested By'},
            {"text":'PI'},
            {"text":'Project Number'},
            {"text": "Gene"},
            {"text":'User Comments'},
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
        srm_list = df_from_ngs_template(template)
        
        #append table to bottom of frame
        #loop through index
        #then loop through element

        for srm in srm_list:
            for entry in srm:
                self.srm_order = srm[0]
                self.requested_by = srm[1]
                self.PI = srm[2]
                self.project_number = srm[3]
                self.gene = srm[4]
                self.user_comments = srm[5]

                
            self.data.append((self.srm_order,
                              self.requested_by,
                              self.PI,
                              self.project_number,
                              self.gene,
                              self.user_comments,
            ))
            
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        def _error_box(self):
            
            self.ngs_date_box.delete(0,END)
            Messagebox.show_warning("Did not find any files. Check SRM's and NGS Date", "I'm sorry Jon.", alert=True)


        def _get_ngs_date(self):
            
            ngs_date = self.ngs_date_box.get()
            
            print(f"NGS Date: {ngs_date}")
            return ngs_date
             
        def _get_subject_line(ngs_date, gene, srm_number):
            
            sub_line = f"NGS {ngs_date} {gene} SRM Order# {srm_number}"
                    
            return sub_line
        
        #pass email object, cage num
        def _get_attachment(email, srm_number, ngs_date):
                        
            attachment_list=[]
            
            ngs_dir = f"Z:\ResearchHome\Groups\millergrp\home\common\\NGS\{ngs_date}\joined"

            #gets all excel files that match the SRM number
            srm_excel_list = glob.glob(ngs_dir+f"/**/*{srm_number}**",recursive=True)
            srm_dir_text_list = []
            
            #gets dir names of each excel file and finds the 1 result text file in the dir
            for excel in srm_excel_list:
                dir_path = os.path.dirname(excel)
                
                #have to choose index of glob.glob.  If no .txt, exception is thrown.  Happens in case of base editing.
                try:
                    srm_dir_text_list.append(glob.glob(dir_path+"/*.txt")[0])
                except:
                    print(f"No text file found in folder: {dir_path}")
                    None

            #can have analysis files but no text files
            if len(srm_excel_list) and len(srm_dir_text_list) != 0:
            #individually attach each file to the attachment list to keep the list flat           
                for file in srm_excel_list:
                    attachment_list.append(file)
                for file in srm_dir_text_list:
                    attachment_list.append(file)
                    
                attachment_list.sort()
                
                for file in attachment_list:
                    email.Attachments.Add(file)
                
                return email
            
            #no excel or text file found. 
            else:
                _error_box(self)
                return None

        def _body_builder(requested_by, gene, project_num):
            
            #pi = pi.split(", ")[1]
            requested_by = requested_by.split(", ")[1]
            
            body=f"""
            <font face="Calibri, Calibri, monospace">
            Hi {requested_by},
            <br><br>
            Attached is the NGS data for your {gene} project.  The CAGE number for this site is {project_num}.
            <br><br>
            The attached .xlsx file(s) include a summary of your data, and the .txt file(s) contain the sequencing reads, as well as details about how the summary sheet was generated. 
            We highly encourage investigators to align their sequencing reads to their gene of interest to verify their results. 
            For more information about our data analysis process and how to interpret the results, check out our video: CRIS.PY Tutorial
            Please let us know if you have any questions.
            <br>
            Thanks
            <br><br>
            </font>               
            """
                
            return body
        
        #gets only rows shown in table and acceses those to generate emails
        intact_rows = self.table.get_rows(visible=True)        
        srm_entries=[]
        
        for row in intact_rows:
            srm_entries.append(row.values)
                
        '''
        'SRM Order #',0
        'PI',1 #052024-Waiting on RIS to add PI column.  all other indicies -= 1
        'Requested By',2
        'Project Number',3
        'Gene', 4
        'User Comments',5
        '''    
        sig = parse_signature()
        
        ngs_date = _get_ngs_date(self)

        #self data is a list of list.  loop through each entry to access each field
        for entry in srm_entries:
            
            srm_number, requested_by, pi, project_num, gene, user_comment = entry
            
            #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            email_recip = [requested_by]
            email_cc = [pi,"Miller, Shondra"]
            
            email_sub = _get_subject_line(ngs_date, gene, srm_number)

            email = _get_attachment(email,srm_number, ngs_date)

            body = _body_builder(requested_by, gene, project_num)

            email.To = ";".join(email_recip) #Requester
            email.CC = ";".join(email_cc).replace(".","") #shondra and PI

            email.Subject = email_sub

            #find html signature file in each individual userprofile
            sig = parse_signature()
            
            email.HTMLBody = body + sig
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)


if __name__ == "__main__":
    import _emailer_gui_RUN_THIS_SCRIPT
    _emailer_gui_RUN_THIS_SCRIPT.app.mainloop()