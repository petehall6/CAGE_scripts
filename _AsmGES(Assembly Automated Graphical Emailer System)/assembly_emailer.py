import ttkbootstrap as tbs
from ttkbootstrap.constants import *
from ttkbootstrap.tableview import Tableview
from emailer_functions import (open_file,
                               df_from_template,
                               parse_signature
                            )
import pandas as pd
import os
import win32com.client
import glob
import re


class Assembly_Tab(tbs.Frame):
    def __init__(self, master_window):
        super().__init__(master_window, padding=(20,20))
        self.pack(fill=BOTH, expand=YES)
        self.button_container = tbs.Frame(self)
        self.button_container.pack(fill=X, expand=YES, pady=(15,10))
        
        self.excel_name = tbs.StringVar(value="")
        
        self.srm_order = tbs.StringVar(value="")
        self.project_number = tbs.StringVar(value="")
        self.PI = tbs.StringVar(value="")
        self.requested_by = tbs.StringVar(value="")
        self.scope = tbs.StringVar(value="")
        self.gene = tbs.StringVar(value="")
        self.species = tbs.StringVar(value="")
        self.vector = tbs.StringVar(value="")
        self.grna1 = tbs.StringVar(value="")
        self.grna2 = tbs.StringVar(value="")
        
        
        self.data = []

        self.create_labels()
        self.create_srm_load_btn()
        self.create_gen_emails_btn()
        self.create_clear_btn()
        self.table = self.create_table()

    def create_labels(self):
    
        self.title_lbl = tbs.Label(
            master = self.button_container,
            text = "Assembly Emailer",
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
            bootstyle=INFO,
            width=25
        )
        
        self.srm_load_btn.grid(column=0,row=1, pady=10)

    def create_gen_emails_btn(self):

        self.gen_emails_btn = tbs.Button(
            master = self.button_container,
            text = "Create Emails",
            command = self.generate_emails,
            bootstyle = WARNING,
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
            {"text":'SRM Order #'},
            {"text":'Project Number'},
            {"text":'PI'},
            {"text":'Requested By'},
            {"text":'Scope'},
            {"text":'Gene'},
            {"text":'Species'},
            {"text":'Vector'},
            {"text":'gRNA1'},
            {"text":'gRNA2'}
        ]

        table = Tableview(
            master = self,
            coldata=columns,
            rowdata=self.data,
            paginated=False,
            searchable=False,
            bootstyle=PRIMARY,
            stripecolor=LIGHT)
        
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
                self.srm_order = srm[0]
                self.project_number = srm[1]
                self.PI = srm[2]
                self.requested_by = srm[3]
                self.scope = srm[4]
                self.gene = srm[5]
                self.species = srm[6]
                self.vector = srm[7]
                self.grna_1 = srm[8]
                self.grna_2 = srm[9]
                    
            self.data.append((
                self.srm_order,
                self.project_number,
                self.PI,
                self.requested_by,
                self.scope,
                self.gene,
                self.species,
                self.vector,
                self.grna_1,
                self.grna_2,
        ))   
        #refresh table with new data.
        self.table.destroy()
        self.table.load_table_data()
        self.table = self.create_table()    

    def generate_emails(self):
        
        signature = parse_signature()
        
        def _email_writer_single(project_details):
            
            srm,project_num, pi, requester, scope, gene, species, vector, grna1, grna2 = project_details
            
            def _get_subject_line():
            
                sub_line = 'Your guides are ready for pickup!'
                return sub_line
            
            def _get_attachment(email, project_num, vector):
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
                
                #find empty vector
                #should pass on "Other"
                try:
                    vector_dir = "Z:\ResearchHome\Groups\millergrp\home\common\Plasmids\Plasmids\DNA Files"
                    vector = re.search(r'\(C.*\)|\(P.*\)', vector)
                    if vector:
                        vector_clean = vector.group().replace("CAGE","").translate(str.maketrans("()","  ")).strip()
                    else:
                        vector_clean = None
                    #*Against my better judgement I am going to hardcode the file names.  I'm sorry future me.  PMH 20250407
                    
                    vector_file_dict = {
                    "SP78":	"SP78-AG65777-BPK1520 w LacZ insert_EMPTY.dna",
                    "SNP28": "SNP28-LentiPuro-AG52963 Empty.dna",
                    "PC292": "PC292_AG96925_pXPR_050_Lenti_empty.dna",
                    "PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
                    "PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
                    "SNP36": "SNP36 LentiCRISPRv2 -AG52961 Empty.dna",
                    "SS363": "SS363 pAW12.lenti.GFP _#104374_empty.dna",
                    "SNP112": "SNP112_AG82416_LentiCRISPRv2GFP Empty.dna",
                    "SP16":	"SP16 - LRG (Lenti_sgRNA_EFS_GFP) - AG65656_empty.dna",
                    "PC529": "PC529_LVA_plent-gRNA-Ametrine_Hongbo_lab_EMPTY_20Ns_update.dna",
                    "PC329": "PC329_pLentiCRIPSRv2_mCherry_AG99154_empty.dna",   
                    "PC226": "PC226-LMA_Hongbo_Empty.dna",
                    "SP76":	"SP76-LV04 U6gRNA PPB Sanger Crispr Vector_MP_newEMPTY_MP.dna",
                    "JK129": "JK129_AG61427_lenti sgRNA-MS2-Zeo_Empty.dna"
                    }
                    
                    vector_path = vector_file_dict.get(vector_clean)
                    
                    vector_file = os.path.join(vector_dir,vector_path)
                    

                except:
                    print("Problem attaching vector file")
                    print(f"Vector choice: {vector}")
                    vector_path = None
                
                
                if latest_ppt is not None:
                    email.Attachments.Add(latest_ppt)
                    
                if vector_path is not None:
                    email.Attachments.Add(vector_file)
                
                return email
        
            def _body_builder(greeting, scope):
                if scope.lower() == "grna - no validation":

                    body=f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Your guides have been cloned and sequence confirmed. They are ready for pick up in M4160. When you enter, find the two mini freezers on the west wall. 
                    Open the right hand freezer and your guides will be in the plastic accordion folder on the door, under the first letter of the PI’s last name.  
                    <br><br>
                    I have included an empty version of the backbone for your reference.
                    <br><br>
                    <br><br>
                    Please let me know if you have any questions.
                    <br><br>
                    Have a good day!
                    <br><br>
                    <br><br>
                    </font>                
                    """
                    
                elif scope.lower() == "grna with validation":
                    body=f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Your guides have been cloned, sequence confirmed, and validated. They are ready for pick up in M4160. When you enter, find the two mini freezers on the west wall.
                    Open the right hand freezer and your guides will be in the plastic accordion folder on the door, under the first letter of the PI's last name.
                    <br><br>
                    I have included an empty version of the backbone for your reference.
                    <br><br>
                    Please let me know if you have any questions.
                    <br><br>
                    Have a good day!
                    <br><br>
                    <br><br>
                    </font>
                    """
                return body
            
        #mail object generator
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            recip_list = [requester]
            email_recip = list(set(recip_list))

            greeting = (f"Hi {requester.split(',')[1]}").replace(",,",",")
                            
            email_sub = _get_subject_line()

            email = _get_attachment(email,project_num,vector)

            body = _body_builder(greeting,scope)
            
            email_cc = [pi,"Miller, Shondra"]

            email.To = ";".join(email_recip)
            email.CC = ";".join(email_cc)

            email.Subject = email_sub

            #find html signature file in each individual userprofile
            
            email.HTMLBody = body + signature
            #Display(False) loads all emails at once and gives focus back to ttk window
            email.Display(False)

        def _email_writer_multi(project_df):
            
            projects = project_df.values.tolist()
            
            #initizlie all the list at once
            
            srm,project_num, pi, requester, scope, gene, species, vectors, grna1, grna2 = ([] for i in range(10))
            
            #project_num, srm, pi, requester, scope, species, batch, gRNA1, gRNA2 = ([] for i in range(9))
            #unpacked nested list into individual list
            srm,project_num, pi, requester, scope, gene, species, vectors, grna1, grna2 = map(list,zip(*projects))
            
            vectors = list(set(vectors))  #remove duplicates from vector list
            
            def _get_attachment(email, project_num, vectors):
                #find powerpoint
                ppt_attachment_list = []
                
                for proj in project_num:
                
                    try:
                        path = "Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS/"
                        for name in glob.glob(os.path.join(path, "*{}".format(proj))):
                            folder = name

                        os.chdir(folder)
                        ppt_list = glob.glob("*.pptx")
                        latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)
                        print(latest_ppt)
                        ppt_attachment_list.append(latest_ppt)
                    
                    except:
                        print("couldn't find slidedeck in CORE Project folder")
                        print("Project Number = {}".format(proj))
                        latest_ppt = None


                        
                for vector in vectors:
                    vector_attachment_list = []
                    try:
                        vector_dir = "Z:\ResearchHome\Groups\millergrp\home\common\Plasmids\Plasmids\DNA Files"
                        vector = re.search(r'\(C.*\)|\(P.*\)', vector)
                        if vector:
                            vector_clean = vector.group().replace("CAGE","").translate(str.maketrans("()","  ")).strip()
                        else:
                            vector_clean = None
                        #*Against my better judgement I am going to hardcode the file names.  I'm sorry future me.  PMH 20250407
                        
                        vector_file_dict = {
                        "SP78":	"SP78-AG65777-BPK1520 w LacZ insert_EMPTY.dna",
                        "SNP28": "SNP28-LentiPuro-AG52963 Empty.dna",
                        "PC292": "PC292_AG96925_pXPR_050_Lenti_empty.dna",
                        "PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
                        "PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
                        "SNP36": "SNP36 LentiCRISPRv2 -AG52961 Empty.dna",
                        "SS363": "SS363 pAW12.lenti.GFP _#104374_empty.dna",
                        "SNP112": "SNP112_AG82416_LentiCRISPRv2GFP Empty.dna",
                        "SP16":	"SP16 - LRG (Lenti_sgRNA_EFS_GFP) - AG65656_empty.dna",
                        "PC529": "PC529_LVA_plent-gRNA-Ametrine_Hongbo_lab_EMPTY_20Ns_update.dna",
                        "PC329": "PC329_pLentiCRIPSRv2_mCherry_AG99154_empty.dna",   
                        "PC226": "PC226-LMA_Hongbo_Empty.dna",
                        "SP76":	"SP76-LV04 U6gRNA PPB Sanger Crispr Vector_MP_newEMPTY_MP.dna",
                        "JK129": "JK129_AG61427_lenti sgRNA-MS2-Zeo_Empty.dna"
                        }
                        
                        vector_path = vector_file_dict.get(vector_clean)
                        
                        print(f"Vector path: {vector_path}")
                        
                        vector_file = os.path.join(vector_dir,vector_path)
                        
                        vector_attachment_list.append(vector_file)
                    

                    except:
                        print("Problem attaching vector file")
                        print(f"Vector choice: {vector}")
                        vector_path = None
                
                print(f"ppt_attachment_list: {ppt_attachment_list}")
                print(f"vector_attachment_list: {vector_attachment_list}")
                
                
                if ppt_attachment_list is not None:
                    for ppt in ppt_attachment_list:
                        email.Attachments.Add(ppt)
                    
                if vector_attachment_list is not None:
                    for vector in vector_attachment_list:
                        email.Attachments.Add(vector)
                    
                return email
        
            def _bullet_maker(srm_order_num):
                bullet_list =""
                
                
                
                for order in zip(list(set(srm_order_num))):
                    order = str(order).replace(",","").translate(str.maketrans("()","  ")).strip()
                    print(f"order {order}")
                    bullet_list += (f"<li>SRM: {order} </li>")
                
                print(f"The bullet_list {bullet_list}")
                    
                return bullet_list
            
            def _body_builder_multi(srm_order_num,greeting, scope):
                
                bullets = _bullet_maker(srm_order_num)
                if scope[0].lower() == "grna - no validation":
                
                    body = f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Your guides have been cloned and sequence confirmed. They are ready for pick up in M4160. When you enter, find the two mini freezers on the west wall. 
                    Open the right hand freezer and your guides will be in the plastic accordion folder on the door, under the first letter of the PI's last name.  
                    <br><br>
                    <ul>
                    {bullets}
                    </ul>
                    <br><br>
                    I have included an empty version of the backbone for your reference.
                    <br><br>
                    Have a good day!
                    <br><br>
                    </font>
                    """
                elif scope[0].lower() == "grna with validation":
                    body = f"""
                    <font face="Calibri, Calibri, monospace">
                    {greeting},
                    <br><br>
                    Your guides have been cloned, sequence confirmed, and validated. They are ready for pick up in M4160. When you enter, find the two mini freezers on the west wall. 
                    Open the right hand freezer and your guides will be in the plastic accordion folder on the door, under the first letter of the PI's last name.
                    <br><br>
                    <ul>
                    {bullets}
                    </ul>
                    <br><br>
                    I have included an empty version of the backbone for your reference.
                    <br><br>
                    Have a good day!
                    <br><br>
                    </font>
                    """

                return body
            
            outlook = win32com.client.Dispatch("Outlook.Application")
            email = outlook.CreateItem(0)
            
            #removes duplicates and rephrases the greeting to a single person
            
            recip_list = requester
            email_recip = list(set(recip_list))
            
            
            print(f"Email receip list: {email_recip}")
            
            if len(email_recip) > 1:
                first_names=[]
                for req in requester:
                    first_names.append(req.split(',')[1])

                first_names = list(set(first_names))
                #insert 'and' in front of last element            
                first_names[-1] = ' and '+first_names[-1] 
     
                greeting = f"Hi {(','.join(first_names))}"
                
            else:
                greeting = f"Hi {email_recip[0].split(',')[1]}"
                            
            email_sub = "Projects ready for pickup"

            email = _get_attachment(email,project_num, vectors)

            body = _body_builder_multi(srm,greeting, scope)

            email.To = ";".join(email_recip)
            
            email_cc = [pi[0],"Miller, Shondra"]
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
                    'Project Number',
                    'PI',
                    'Requested By',
                    'Project Scope',
                    'Gene',
                    'Species',
                    'Vector',
                    'gRNA1',
                    'gRNA2'
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