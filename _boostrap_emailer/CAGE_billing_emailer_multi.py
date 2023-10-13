# -*- coding: utf-8 -*-
"""

07-24-2020
@author: Rachel Levine
Program to create emails for letting users know they can pick up their complete cell line project
Requires input of excel export from SRM2 from Billing page of CAGE PROJECTS
python 2.7
runs with ctrl+F5

"""


import pandas as pd
import numpy as np
import os
import sys
import glob
from datetime import date


os.chdir(os.path.dirname(os.path.abspath(__file__)))


def create_bill_df():
    try:
        list_of_files = glob.glob("CAGEServices_Excel Export*.xls")
        latest_file = max(list_of_files, key=os.path.getctime)
        print(latest_file)
    except:
        print(
            "No export file found.  Make sure this python file is saved in the same folder with your excel export"
        )
        sys.exit()
    bill_df = pd.read_excel(
        latest_file, sheet_name=0
    )  # opens first sheet of excel file
    bill_df.to_csv("email.csv", encoding="utf-8-sig")
    bill_df = pd.read_csv("email.csv")

    bill_df.rename(
        columns={
            "SRM Project #": "SRM_no",
            "Requested By": "user",
            "Species": "species",
            "Project Number": "proj_no",
            "Target Gene Name": "gene",
            "Project Objective": "objective",
            "Cell Line of Choice": "cell_line",
            "Project Scope": "scope",
        },
        inplace=True,
    )

    bill_df = bill_df.reset_index(drop=True)
    return bill_df


def text_fxn(user, PI, cell_line, gene, objective, scope):
    PI_name = PI.split(", ")[1]
    user_name = user.split(", ")[1]
    if scope == "Edited cell pool":
        text = (
            "Hi "
            + PI_name
            + " and "
            + user_name
            + ",<br><br>Great news! Your "
            + gene
            + "_"
            + cell_line
            + " edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.<br><br>The last slide is the most informative.  We were able to get over XX% total editing in the pool with ~XX% out of frame indels.<br><br>We have a contactless pickup system in place right now.  Please coordinate with XXXXX to let XXXX know a good time window for you to pick up these cells. During the agreed upon time, XXXX will place the frozen vials of cells in a dry ice bucket in M4170. The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  If you would like the live cultures as well, please come in the next day or so.  The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.<br><br>Don't hesitate to contact me if you have any questions.<br>Best,<br><br>SM"
        )
    elif scope == "Cell Fitness/Dependency Assay":
        text = f"Hi {PI_name} and {user_name},<br><br>Great news! Your {gene} {cell_line} fitness assay is complete. Please see the attached slide deck for details.<br><br>We do/do not see a strong dependency for this gene in this background.<br><br>Please let me know if you have any questions.<br><br>Best,<br>SM"
    else:
        text = (
            "Hi "
            + PI_name
            + " and "
            + user_name
            + ",<br><br>Great news! Your "
            + cell_line
            + "_"
            + gene
            + "_"
            + objective
            + " project is complete and ready for pick up.  Please see the attached slide deck for details.<br> We currently have a contactless pickup system in place.  Please arrange a time window with XXXXX in which someone can pick up the cells.  At the agreed upon time, he/she will place your frozen vials of cells into a dry ice bucket in M4170.  The dry ice bucket will be on the counter in front of you as you walk in.  Your live cultures will be in the first incubator to the right (top incubator, bottom shelf) and labeled accordingly. Please also bring dry ice for the pickup.<br><br>As always, please let me know if you have any questions.<br><br>Best,<br>SM"
        )
    return text


def fetch_attachment(proj_no):
    try:
        path = "Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS/"
        for name in glob.glob(os.path.join(path, "*{}".format(proj_no))):
            folder = name
        # print("CORE Project folder is called {}".format(os.path.basename(folder)))
        os.chdir(folder)
        ppt_list = glob.glob("*.pptx")
        latest_ppt = folder + "/" + max(ppt_list, key=os.path.getctime)

        # print("Slidedeck path = {}".format(path1))
        return latest_ppt
    except:
        print("couldn't find slidedeck in CORE Project folder")
        print("Project Number = {}".format(proj_no))
        return None


def add_attachments(mail_obj, proj_list):
    for proj_no in proj_list:
        attachment = fetch_attachment(proj_no)
        if attachment is not None:
            mail_obj.Attachments.add(attachment)

    return mail_obj


def Emailer(text, signature, subject, recipient, recipient2, proj_list, multi=False):
    import win32com.client as win32

    outlook = win32.Dispatch("outlook.application")
    mail = outlook.CreateItem(0)
    mail.To = recipient + ";" + recipient2
    # mail.cc = recipient2
    mail.bcc = "Porter, Shaina"
    mail.Subject = subject
    mail.HtmlBody = text + signature
    mail.Display(False)
    mail = add_attachments(mail, proj_list)


def gen_pi_df(pi_name, bill_df):
    pi_df = bill_df[bill_df["PI"] == pi_name]
    # proj_list = pi_df['proj_no'].unique().tolist()
    return pi_df


def handle_pi_df(pi_df):
    if len(pi_df) > 1:
        email_multi_project(pi_df)
        return
    else:
        email_single_project(pi_df)
        return


def email_single_project(pi_df):
    return


def email_multi_project(pi_df):
    def group_by_scope():
        return

    def bullet_point_text(row):
        return_text = f"<li>SRM {row['SRM_no']}: {row['gene']} {row['objective']} in {row['cell_line']} ({row['scope']})</li>"
        return return_text

    list_of_bullet_points = pi_df.apply(bullet_point_text, axis=1)
    return list_of_bullet_points


def original_emailer(bill_df, signature):
    for i in range(bill_df["SRM_no"].count()):
        SRM_no = bill_df.loc[i, "SRM_no"]
        PI = bill_df.loc[i, "PI"]
        user = bill_df.loc[i, "user"]
        gene = bill_df.loc[i, "gene"]
        cell_line = str(bill_df.loc[i, "cell_line"])
        objective = bill_df.loc[i, "objective"]
        scope = bill_df.loc[i, "scope"]
        try:
            proj_no = bill_df.loc[i, "proj_no"].split(" ")[0]
        except:
            print(
                "Please check that your project numbers are entered correctly in the CAGEServices_Excel Export"
            )
            sys.exit()
        # print(cell_line+gene+objective)
        text = text_fxn(user, PI, cell_line, gene, objective, scope)
        if scope == "Edited cell pool":
            subject = gene + " " + cell_line + " edited cell pool complete"
        elif scope == "Cell Line Creation":
            subject = cell_line + " " + gene + " " + objective + " cell line complete "
        elif scope == "Cell Fitness/Dependency Assay":
            subject = f"Cell Fitness/Dependency Assay for {gene} in {cell_line}"
        Emailer(str(text), signature, subject, user, PI, proj_no)

    print("All draft emails are created")
    return


def main():
    signature = "<br><br>Shondra Miller, Ph.D. <br> Director, Center for Advanced Genome Engineering (CAGE)<br>Associate Director, Shared Resources for Comprehensive Cancer Center<br>Associate Member, Department of Cell and Molecular Biology<br>St. Jude Children's Research Hospital<br>262 Danny Thomas Place <br>Memphis, TN  38105<br>Office: (901)595-7313 - Rm M4413<br>Lab: (901)-595-7314 - Rm M4160<br>"

    bill_df = create_bill_df()
    original_emailer(bill_df, signature)
    return


if __name__ == "__main__":
    main()
