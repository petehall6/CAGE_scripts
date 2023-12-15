def _body_builder(pi,requested_by,sample_type,success,edit,gene,cage_number):
        
        #TODO work out greeting
        greeting = f"""Hello, {pi.split(",")[1]} and {requested_by.split(",")[1]}"""
        
        if sample_type =="Tail Snip/Toe Snip":
            if success.upper() == "YES":
                if edit.upper() == "CKO":
                    print("CKO Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} {edit} animals are positive for both the 5' and 3' loxP sites.  
                    The CAGE numbers for these sites are {cage_number}.
                    <br><br>
                    I have also included the data for large deletions between the two guide sites.  Deletion animals could be used to generate a 
                    germline KO if bred to homozygosity and viable.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """
                elif edit.upper() == "KO":
                    print("KO Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} {edit} animals contain out of frame indels (highlighted in green).
                    The CAGE numbers for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """
                elif edit.upper() == "KI":
                    print("KI Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} {edit} animals contain both the 5' and 3' junctions between your desired integration and your target site..
                    The CAGE number for this site is {cage_number}.  I've also attached the data for the WT site."
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """
                elif edit.upper() == "DEL":
                    print("DEL Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} deletion animals contain deletions.
                    The CAGE numbers for this site are {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """   
                elif edit.upper() == "SSODN":
                    print("SSODN \n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} ssODN animals contain your desired mutation.
                    The CAGE number for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """ 
                elif edit.upper() == "POINT MUTATION":
                    print("POINT mutation \n\n")
                    body=f"""{greeting},
                    <br><br>
                    Great news! {submitted_num} of the {gene} point mutation animals contain your desired mutation.
                    The CAGE number for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """
                else:
                    print("YES ELSE\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} project.
                    The CAGE number for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """ 
                return body
                
            elif success.upper() == "NO":
                if edit.upper() == "CKO":
                    print("NO - CKO Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} {edit} project.  Unfortunately, none of these animals contain the desired loxP sites.
                    The CAGE numbers for these sites are {cage_number} 5' and 3'.
                    I have also included the data for large deletions between the two guide sites. These animals could be used to generate a germline KO if bred to homozygosity and viable.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """
                elif edit.upper() == "KO":
                    print("NO - KO Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} {edit} project.  Unfortunately, none of these animals show any editing.
                    The CAGE numbers for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>               
                    """
                elif edit.upper() == "KI":
                    print("NO - KI Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} {edit} project.  Unfortunately, none of the animals contain both the 5' and 3' junctions 
                    between your desired integration and your target site.
                    The CAGE numbers for this site is {cage_number}.  I've also attached the data for the WT site.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>               
                    """
                elif edit.upper() == "DEL":
                    print("NO - DEL Project\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} {edit} project.  Unfortunately, none of tese animals contain deletions.  
                    The CAGE numbers for this site are {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """   
                else:
                    print("NO - ELSE\n\n")
                    body=f"""{greeting},
                    <br><br>
                    Attached is the NGS data for your {gene} {edit} project.
                    The CAGE number for this site is {cage_number}.
                    <br><br>
                    Don't hesitate to contact me if you have any questions.
                    <br><br>
                    Best,
                    <br><br>
                    <br><br>                
                    """ 
                return body
        
        else:
            body=f"""{greeting},
                        <br><br>
                        Attached is the NGS data for your {gene} project.
                        The CAGE number for this site is {cage_number}.
                        <br><br>
                        Don't hesitate to contact me if you have any questions.
                        <br><br>
                        Best,
                        <br><br>
                        <br><br>
                        """
            return body

def _email_writer(pi, requested_by, department, gene, edit, edit_size, injection_core, cage_number, ngs_date, success_num, submitted_num, notes, srm_number, success):
            
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
                  #  print("couldn't find slidedeck in CORE Project folder")
                  #  print("Project Number = {}".format(project_num))
                    latest_ppt = None

                if latest_ppt is not None:
                    email.Attachments.Add(latest_ppt)
                    
                
                
                return email
        
            def _body_builder(edit, success, greeting, scope, cell_line, objective, line_lead):
                if edit.upper() == "CKO":

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
            
            signature = parse_signature()
            
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
    
table_rows = [('2665904', 'Kanneganti, Thirumala-Devi', 'Baskaran, Yogi', 'Baskaran, Yogi', 'CAGE84', 'RIPK3-mBFP2_KI', '18', 'Well of Plate', 'Tail Snip/Toe Snip', 'CBT', 'NO', '1', '4', 'KO', '5', 'TCU', '111423', ['CAGE84_hDGCR6L_F_R_short'],'testing_notes')]

for row in table_rows:
    proj_data = list(row)
    
    print(proj_data)
    #unpack proj_data
    (srm_number,
    pi,
    requested_by,
    entered_by,
    cage_number,
    gene,
    sample_num,
    sample_format,
    sample_type,
    department,
    success,
    success_num,
    submitted_num,
    edit,
    edit_size,
    injection_core,
    ngs_date,
    cage_programs,
    notes) = proj_data
    
    #update excel sheet
body = _body_builder(pi,requested_by,sample_type,success,edit,gene,cage_number)
#_email_writer(pi, requested_by, department, gene, edit, edit_size, injection_core, cage_number, ngs_date, success_num, submitted_num, notes, srm_number, success)

print(body)

