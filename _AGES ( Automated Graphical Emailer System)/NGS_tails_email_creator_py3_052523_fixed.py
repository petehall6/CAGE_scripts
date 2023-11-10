# -*- coding: utf-8 -*-

""" author: Rachel Levine 05/2020
run with python 3.7
run with ctrl+ F5 or "run python file from terminal"  in VS Code """


import pandas as pd
#import numpy as np
from datetime import date
import os
import glob
import sys
import openpyxl
import getpass
import win32com.client as win32
import math


username = getpass.getuser()

# file_date =  str(date.today().strftime('%m_%d_%Y'))
os.chdir(os.path.dirname(os.path.realpath(__file__)))
curr_dir = os.getcwd()
print(curr_dir)

now = date.today()

column_names = [
    "PI",
    "dept",
    "gene",
    "mod",
    "size",
    "Start date",
    "date",
    "Days to completion",
    "injection core",
    "Notes",
    "proj_no",
]

report_df = pd.DataFrame(columns = column_names)

# department abbreviations
dept_dict = {
    "CHEMICAL BIOLOGY AND THERAPEUTIC": "CBT",
    "CELL & MOLECULAR BIOLOGY": "CMB",
    "COMPUTATIONAL BIOLOGY": "Comp Bio",
    "DEVELOPMENTAL NEUROBIOLOGY": "DNB",
    "PHARMACEUTICAL SCIENCES": "Pharm Sciences",
    "STRUCTURAL BIOLOGY": "Structural Bio",
    "TUMOR BIOLOGY": "Tumor Bio",
    "BONE MARROW TRANS & CELL THRPY": "BMT",
    "INFECTIOUS DISEASES": "Infections Dis",
}

# find the latest excel export from SRM

list_of_files = glob.glob(f"C:\Users\jklein1\Downloads\TargetedDeepSequencing_Excel Export*.xls")

#For debugging
#list_of_files = glob.glob("TargetedDeepSequencing_Excel Export*.xls")

try:
    latest_file = max(list_of_files, key=os.path.getctime)
except:
    print(
        "No export file found.  Make sure this file is saved in the same folder with your SRM excel export file."
    )
    sys.exit()
# Ask user for the NGS date where all the data is currently saved
NGS_date = input("Enter the NGS date ").strip()
user = input(
    "Enter your initials "
).strip()  # ask the user who they are so the program knows whose signature to include
if user == "RL":
    signature = "<br><br>Rachel M Levine, Ph.D.<br>Center for Advanced Genome Engineering (CAGE)<br>St. Jude Children's Research Hospital<br>Department of Cell and Molecular Biology<br>262 Danny Thomas Place<br>Memphis, TN  38105<br>Office: (901)595-7314 - Rm M4160"
elif user == "SNP":
    signature = "<br><br>Shaina Porter, Ph.D.<br>Center for Advanced Genome Engineering (CAGE)<br>St. Jude Children's Research Hospital<br>Department of Cell and Molecular Biology<br>262 Danny Thomas Place<br>Memphis, TN  38105<br>Office: (901)595-7314 - Rm M4160"
elif user == "JK":
    signature = "<br><br>Jonathon Klein<br>Center for Advanced Genome Engineering (CAGE)<br>St. Jude Children's Research Hospital<br>Department of Cell and Molecular Biology<br>262 Danny Thomas Place<br>Memphis, TN  38105<br>Office: (901)595-7314 - Rm M4160"
else:
    signature = input(
        "You must be new, please paste in what you would like as your email signature:"
    ).strip()

NGS_df = pd.read_excel(latest_file, sheet_name=0)  # opens first sheet of excel file
NGS_df.to_csv("NGS.csv", encoding="utf-8-sig")
NGS_df = pd.read_csv("NGS.csv")
# this block updates the Completed animal models file if this project is a mouse model that is successful
def update_mice(report_df):
    # open existing mouse models records
    print("update_mice")
    projects_df = pd.read_excel(
        rf"C:\Users\{username}\St. Jude Children's Research Hospital\Team-CAGE - General\Completed Animal Models (CAGE).xlsx",
        sheet_name=0,
    )
    projects_df.to_csv("mice.csv", encoding="utf-8-sig")
    projects_df = pd.read_csv("mice.csv")
    projects_df.rename(
        columns={
            "Investigator": "PI",
            "Department": "dept",
            "Gene name": "gene",
            "Project type": "mod",
            "insert size": "size",
            "End date": "date",
            "Initials of who completed project": "injection core",
            "Principal Investigator": "PI",
            "CAGE Project #": "proj_no",
        },
        inplace=True,
    )
    projects_df = projects_df.reset_index(drop=True)
    # update mouse models records
    projects_df = projects_df.append(report_df, sort=False, ignore_index=True)
    wb = openpyxl.load_workbook(
        f"C:\Users\{username}\St. Jude Children's Research Hospital\Team-CAGE - General\Completed Animal Models (CAGE).xlsx"
    )
    ws = wb.worksheets[0]
    # adds the project info from report_df to the next line in completed animal models
    try:
        ws = ws.append(report_df)
        print("added line to ws")
    except:
        print("No data to add to the Completed CAGE animal projects list")
    wb.save(
        f"C:\Users\{username}\St. Jude Children's Research Hospital\Team-CAGE - General\Completed Animal Models (CAGE).xlsx"
    )

# renames the colums in the df created from the excel export to easier to use variables
NGS_df.rename(
    columns={
        "CAGE Project #": "proj_no",
        "Gene Name/Gene ID": "gene",
        "Number of Tube Samples/Plates": "no_sample",
        "Sample Format": "format",
        "Sample Type": "s_type",
        "SRM Order #": "order",
        "Requested By": "requester",
        "Principal Investigator": "PI",
        "Entered By": "enter"
        # "Number of positive animals": "number"
    },
    inplace=True,
)

NGS_df = NGS_df.reset_index(drop=True)

save_df = pd.DataFrame()
skip_df = pd.DataFrame()
# this block finds all the folders corresponding to the project number in the indeicated NGS date, and allows the user to choose in which folders to find the attachments
def files(proj_no, edit):
    dir_list = []
    attach_list = []
    latest_xlsx = []
    latest_txt = []
    print("Here are the folders containing " + proj_no)
    dir_list = glob.glob(
        "Z:\ResearchHome\Groups\millergrp\home\common\\NGS\\"
        + NGS_date
        + "\\joined\\"
        + proj_no
        + "*"
        + "\\"
    )
    for path in dir_list:
        print(str(dir_list.index(path) + 1) + " " + path)
    while True:
        try:
            dir_choice = input(
                "Enter the number of the directories where your files are saved.  If more than 1, separate each by a space (i.e. 1 3 5)"
            ).strip()
            for ans in dir_choice.split():
                # print(dir_list[int(ans)-1])
                attach_list.append(dir_list[int(ans) - 1])
        except IndexError:
            print("invalid choice, try again")
            attach_list = []
        else:
            break
    print(attach_list)
    for folder in attach_list:
        os.chdir(folder)
        list_of_files = glob.glob("*.xlsx")
        try:
            latest_xlsx.append(folder + max(list_of_files, key=os.path.getctime))
            list_of_files = glob.glob("*.txt")
            latest_txt.append(folder + max(list_of_files, key=os.path.getctime))
        except:
            print("no files found in folder " + folder)
    return latest_xlsx, latest_txt


# this block composes the email, using the email body text, the subject, recipients, and paths for the attachments
def Emailer(xlsx, txt, text, subject, recipient, recipient2, geno_advice):


    outlook = win32.Dispatch("Outlook.Application")
    mail = outlook.CreateItem(0)
    mail.To = recipient
    mail.cc = "Miller, Shondra;" + recipient2
    mail.Subject = subject
    mail.HtmlBody = text
    mail.Display(False)
    # this is where the files are attached.  It is a finicky command
    
    print(geno_advice)
    
    
    if geno_advice != "Nan":
    
        try:
            mail.Attachments.Add(geno_advice)
        except:
            print("no genotyping advice to add")
            
            
    for ind in range(len(xlsx)):
        
        try:
            attachment1 = xlsx[ind]
            attachment1 = str(attachment1).replace("\\", "\\\\")
            mail.Attachments.Add(attachment1)
        except Exception as err:
            print(err)
            print("couldn't find attachments")
        try:
            attachment2 = txt[ind]
            mail.Attachments.Add(attachment2)
        except Exception as err:
            print(err)
            print("couldn't find attachments")
    return None


# this block creates the correct text for the email body, based on the project type and whether or not the project was successful
def text_fxn(
    gene, edit, greeting, results, proj_no, s_type, number, advice_msg
):  # number, gene, edit, PI, results, proj_no):
    Crispy = "<br><br>The attached .xlsx file(s) include a summary of your data, and the .txt file(s) contain the sequencing reads, as well as details about how the summary sheet was generated.  We highly encourage investigators to align their sequencing reads to their gene of interest to verify their results.  For more information about our data analysis process and how to interpret the results, check out our video: <a href='https://s.stjude.org/video/player.html?videoId=6000021936001'>CRIS.PY Tutorial</a>"
    closing = "<br>Please let us know if you have any questions.<br>Thanks!"
    if s_type == "Tail Snip/Toe Snip":

        if results == "yes" or results == "y":
            if edit == "CKO":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals are positive for both the 5' and 3' loxP sites. The CAGE numbers for these sites are "
                    + proj_no
                    + " .<br> I have also included the data for large deletions between the two guide sites. Deletion animals could be used to generate a germline KO if bred to homozygosity and viable."
                )
            elif edit == "KO":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals contain out of frame indels (highlighted in green). The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            elif edit == "KI":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals contain both the 5' and 3' junctions between your desired integration and your target site. The CAGE number for this site is "
                    + proj_no
                    + ". I've also attached the data for the WT site."
                )
            elif edit == "deletion":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals contain deletions. The CAGE numbers for this site are "
                    + proj_no
                    + "."
                )
            elif edit == "ssODN":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals contain your desired mutation. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            elif edit == "PM":
                text = (
                    ", <br> Great news! "
                    + number
                    + " of the "
                    + gene
                    + " "
                    + edit
                    + " animals contain your desired mutation. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            elif edit == "data":
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " project. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            else:
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. The CAGE number for this site is "
                    + proj_no
                    + "."
                )

        elif results == "no":
            if edit == "CKO":
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. Unfortunately, none of these animals contain the desired loxP sites. The CAGE numbers for these sites are "
                    + proj_no
                    + " 5' and 3'.<br> I have also included the data for large deletions between the two guide sites. These animals could be used to generate a germline KO if bred to homozygosity and viable."
                )
            elif edit == "KO":
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. Unfortunately, none of these animals show any editing. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            elif edit == "KI":
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. Unfortunately, none of the animals contain both the 5' and 3' junctions between your desired integration and your target site. The CAGE number for this site is "
                    + proj_no
                    + ". I've also attached the data for the WT site."
                )
            elif edit == "deletion":
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. Unfortunately, none of these animals contain deletions. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
            else:
                text = (
                    ", <br> Attached is the NGS data for your "
                    + gene
                    + " "
                    + edit
                    + " project. The CAGE number for this site is "
                    + proj_no
                    + "."
                )
        else:
            text = (
                ", <br> Attached is the NGS data for your "
                + gene
                + " "
                + edit
                + " project. The CAGE number for this site is "
                + proj_no
                + "."
            )
    else:
        text = (
            ", <br> Attached is the NGS data for your "
            + gene
            + " "
            + edit
            + " project. The CAGE number for this site is "
            + proj_no
            + "."
        )
    full_text = greeting + text + Crispy + advice_msg + closing
    return full_text


# this block iterates through all of the projects in the excel export df and displays them to the user to prep to create the text, find the attachments, and compose the email
def send_project(NGS_df):
    column_names = [
        "PI",
        "dept",
        "gene",
        "mod",
        "size",
        "Start date",
        "date",
        "Days to completion",
        "injection core",
        "Notes",
        "proj_no",
    ]
    report_df = pd.DataFrame(columns=column_names)
    geno_advice = "Nan"
    advice_msg = ""
    proj_no = str(NGS_df.loc[i, "proj_no"])
    gene = str(NGS_df.loc[i, "gene"])
    s_type = str(NGS_df.loc[i, "s_type"])
    requester = str(NGS_df.loc[i, "requester"])
    print("Requestor: " + requester)
    order = str(int(NGS_df.loc[i, "order"]))
    PI = str(NGS_df.loc[i, "PI"])

    edit = input(
        "What type of edit is it? Your choices are data (for a generic message), KO, KI, CKO, deletion, ssODN, PM "
    ).strip()
    results = input("Was their project successful? ").strip()
    mod = edit
    if edit == "KO":
        mod = "knockout"
    elif edit == "KI":
        mod = "knockin"
    elif edit == "PM":
        mod = "point mutation"
    elif edit == "CKO":
        mod = "conditional knockout"

    if s_type == "Tail Snip/Toe Snip":
        number = input("How many animals? ").strip()
    else:
        number = 0
    f_name = PI.split(", ")[1]
    l_name = PI.split(", ")[0]
    f_name_req = requester.split(", ")[1]
    # this step automatically adds some extra recipients for projects that are known to include more than one person
    if requester == "Dillard Stroud, Miriam E":
        requester = "Dillard Stroud, Miriam E; Ansari, Shariq"
    if requester == "Zhang, Tina":
        requester = "Dillard Stroud, Miriam E; Zhang, Tina"
    if PI == "Kanneganti, Thirumala-Devi":
        PI = "Kanneganti, Thirumala-Devi; malireddi.subbarao@stjude.org; Nadella, Vinod; Baskaran, Yogi; Chadchan, Sangappa; Sharma, Bhesh; "
    elif PI == "Geiger, Terrence L":
        PI = "Geiger, Terrence L; Alli, Rajshekhar"
    elif PI == "Klco, Jeffery":
        PI = "Klco, Jeffery; Westover, Tamara"
    elif PI == "Kundu, Mondira":
        PI = "Kundu, Mondira; Li-Harms, Xiujie"
    elif PI == "Schuetz, John":
        PI = "Schuetz, John; Wang, Yao"
    elif PI == "Crispino, John":
        PI = "Crispino, John; Hall, Trent"
    elif PI == "Downing, James":
        PI = "Parganas, Evan"
    elif PI == "Hatley, Mark":
        PI = "Hatley, Mark; Garcia, Matthew"
    elif PI == "Chi, Hongbo":
        PI = "Chi, Hongbo; Rankin, Sherri"
    elif PI == "Thomas, Paul":
        PI = "Thomas, Paul; Sisti, Resha; Van De Velde, Lee Ann"
    elif PI == "Yu, Jiyang":
        PI = "Yu, Jiyang; Yang, Xu"
    # this block begins the update mice function and adds the genotyping advice if the submission is from Hartmut or Valerie and it was successful.
    if enter == "Berns, Hartmut" or enter == "Covington, Arthur":
        recipient1 = PI
        recipient2 = "Sublett, Jack; Li, Ling"
        greeting = "Hi Dr. " + l_name
        if results == "yes" or results == "y":
            geno_advice = r"Z:/ResearchHome/Groups/millergrp/home/common/Protocols and SOPs/NGS/tails/CAGE Genotyping Advice.pdf"
            advice_msg = "<br>Additionally, I’m attaching a pdf file that contains some helpful advice for following up on the F0 genotyping results sent today."
            dept = input("What department is this PI in?").strip()
            # replace the full name for the department with the abbrev for the dept
            for key, value in dept_dict.items():
                if key in dept.upper():
                    dept = dept.upper().replace(key, value)
                else:
                    pass
            print(dept)
            size = input(
                "If this is a KI, how many bps is the insert? (NA if it's not)"
            ).strip()
            if size == "":
                size = "NA"
            total = input("How many samples have been submitted so far?").strip()
            add = input("Any notes to add?").strip()
            note = (
                str(number)
                + " positive animals out of "
                + total
                + " in this NGS round "
                + add
            )
            report_df = [
                PI,
                dept,
                gene,
                mod,
                size,
                "",
                now.strftime("%m/%d/%Y"),
                "",
                "HB",
                note,
                proj_no,
            ]
            update_mice(report_df)

    elif enter == "Stewart, Valerie":
        recipient1 = PI
        recipient2 = enter
        greeting = "Hi Dr. " + l_name

        if results == "yes" or results == "y":
            geno_advice = "Z:\\ResearchHome\\Groups\\millergrp\\home\\common\\Protocols and SOPs\\NGS\\Tails\\CAGE Genotyping Advice.pdf"
            advice_msg = "<br>Additionally, I’m attaching a pdf file that contains some helpful advice for following up on the F0 genotyping results sent today."
            dept = input("What department is this PI in?").strip()
            # replace the full name for the department with the abbrev for the dept
            for key, value in dept_dict.items():
                if key in dept.upper():
                    dept = dept.upper().replace(key, value)
                else:
                    pass
            print(dept)
            size = input(
                "If this is a KI, how many bps is the insert? (NA if it's not)"
            ).strip()
            if size == "":
                size = "NA"
            add = input("Any notes to add?").strip()
            total = input("How many samples have been submitted so far?").strip()
            note = (
                str(number)
                + " positive animals out of "
                + total
                + " in this NGS round "
                + add
            )
            report_df = [
                PI,
                dept,
                gene,
                mod,
                size,
                "",
                now.strftime("%m/%d/%Y"),
                "",
                "VS",
                note,
                proj_no,
            ]
            update_mice(report_df)
    else:
        recipient1 = requester
        recipient2 = PI
        greeting = "Hi " + f_name_req
        advice_msg = ""
    file_result = files(proj_no.split(" ")[0], edit)
    xlsx = file_result[0]
    txt = file_result[1]
    if proj_no.split(" ")[0] != proj_no.split(" ")[-1]:
        file_result_2 = files(proj_no.split(" ")[-1], edit)
        xlsx_2 = file_result_2[0]
        txt_2 = file_result_2[1]
        xlsx.extend(xlsx_2)
        txt.extend(txt_2)
    text = text_fxn(
        gene, edit, greeting, results, proj_no, s_type, str(number), advice_msg
    )
    print(type(xlsx))
    print(xlsx)
    print(type(txt))
    print(txt)
    Emailer(
        xlsx,
        txt,
        str(text) + signature,
        "NGS " + NGS_date + " " + gene + " SRM order " + order,
        recipient1,
        recipient2,
        geno_advice,
    )


for i in range(NGS_df["proj_no"].count()):
    proj_no = str(NGS_df.loc[i, "proj_no"])
    print("Project Number: " + proj_no)
    gene = str(NGS_df.loc[i, "gene"])
    s_type = str(NGS_df.loc[i, "s_type"])
    requester = str(NGS_df.loc[i, "requester"])
    print("Requestor: " + requester)
    order = str(int(NGS_df.loc[i, "order"]))
    PI = str(NGS_df.loc[i, "PI"])
    enter = str(NGS_df.loc[i, "enter"])
    print("Entered by: " + enter)
    print("PI: " + PI)
    rows = NGS_df.loc[i, :]
    while True:
        list_df = input(
            "Enter y to send this project, s to save it for later and k to skip "
        ).strip()
        if list_df not in ("y", "s", "k"):
            list_df = input(
                "Enter y to send this project, s to save it for later and k to skip "
            ).strip()
        else:
            break
    if list_df == "y":
        report = send_project(NGS_df)

    elif list_df == "s":
        save_df = save_df.append(rows, ignore_index=True)
        NGS_df.drop([i])
        print(save_df)
    elif list_df == "k":
        skip_df = skip_df.append(rows, ignore_index=True)
        NGS_df.drop([i])
        print(skip_df)
while len(save_df.index) > 0:
    proj_no = save_df.loc[i, "proj_no"]
    print("Project Number: " + proj_no)
    gene = save_df.loc[i, "gene"]
    s_type = save_df.loc[i, "s_type"]
    requester = save_df.loc[i, "requester"]
    print("Requestor: " + requester)
    order = str(int(save_df.loc[i, "order"]))
    PI = save_df.loc[i, "PI"]
    print("PI: " + PI)
    rows = save_df.loc[i, :]
    while True:
        list_df = input("Enter y to send this project, or k to skip ")
        if list_df not in ("y", "k"):
            list_df = input("Enter y to send this project, or k to skip ")
        else:
            break
    if list_df == "y":
        send_project(save_df)
        save_df.drop([i])
    elif list_df == "k":
        skip_df = skip_df.append(rows, ignore_index=True)
        save_df.drop([i])
        print(skip_df)

print("All draft emails are created")

""" if __name__ == "__main__":
    main() """
