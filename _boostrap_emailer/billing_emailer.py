import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import *
import pandas as pd



class Billing_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        #self.colors = master_window.style.colors
        self.excel_name = tbs.StringVar(value="")
        self.srm_order = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.project_scope = tbs.StringVar(value="")
        self.cell_line = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.data = []
        
        self.create_buttonbox()
        self.create_labels()
        
        self.table = self.create_table()


    def create_buttonbox(self):
        button_container = tbs.Frame(self)
        button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.srm_load_btn = tbs.Button(
            master = button_container,
            text = "Select SRM Template",
            command = self.load_srm,
            bootstyle=SUCCESS,
            width=15
        )
        
        self.srm_load_btn.pack(side=LEFT, padx=5)
        
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
            {"text":'Project Scope'},
            {"text":'Cell Line of Choice'},
            {"text":'Project Objective'},
            {"text":'Target Gene Name'}
        ]

        table= Tableview(
            master = self,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=True,
            bootstyle=PRIMARY,
            stripecolor=LIGHT   
        )

        
        table.pack(fill=BOTH, expand=YES, padx=10, pady=10)
        return table

    #gets name of SRM template and sets the excel_name variable
    def load_srm(self):
        #get name of .xls
        template = open_file()
        self.excel_lbl.config(text=template)
        #convert to dataframe
        
        
        #this returns a list of lists.  Will need to figure out how to unpack that for however long the list is
        a = df_from_template(template)
        
        print(f"This is the srm_order column: {a}")
        
        
        print(self.data)
        #append table to bottom of frame



