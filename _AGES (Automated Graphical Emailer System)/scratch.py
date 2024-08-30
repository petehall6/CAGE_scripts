import os
import glob
import win32com.client

def _get_attachment(email, srm_number, ngs_date):
                        
            found_file = False
            attachment_list=[]
            try:
                ngs_dir = f"Z:\ResearchHome\Groups\millergrp\home\common\\NGS\{ngs_date}\joined"

                print(os.getcwd())

                #*Change to for loop.  Currently just finds the first index of glob.glob list
                
                #gets all excel files that match the SRM number
                srm_excel_list = glob.glob(ngs_dir+f"/**/*{srm_number}**",recursive=True)
                srm_dir_text_list = []
                
                input(srm_excel_list)
                
                
                
                #gets dir names of each excel file and finds the 1 result text file in the dir
                for excel in srm_excel_list:
                    dir_path = os.path.dirname(excel)
                    srm_dir_text_list.append(glob.glob(dir_path+"/*.txt"))

                input(srm_dir_text_list)
                found_file = True
                

                
                
                #individually attach each file to the attachment list to keep the list flat
                for file in srm_excel_list:
                    attachment_list.append(file)
                #for file in srm_dir_text_list:
                   # attachment_list.append(file)

            except:
                print("couldn't find excel file or text file.  Check SRM#'s")
                found_files = False

            if found_file == True:
                for file in attachment_list:
                    email.Attachments.Add(file)
                
            
            return email

        
outlook = win32com.client.Dispatch("Outlook.Application")

email = outlook.CreateItem(0)
srm_number = 829208
ngs_date = "061824"



        
email = _get_attachment(email, srm_number, ngs_date)

email.Display(False)