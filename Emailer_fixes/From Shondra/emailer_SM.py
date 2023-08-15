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
import win32com.client as win32


CORE_PROJECTS_DIR = r"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS"


def gen_bill_df():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
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
    bill_df.to_csv("email.csv", index=False)
    bill_df = pd.read_csv("email.csv")
    bill_df.rename(
        columns={
            "SRM Project #": "SRM_no",
            "Requested By": "user",
            "Species": "species",
            "Project Number": "proj_no",
            "Target Gene Name": "gene",
        },
        inplace=True,
    )
    try:
        bill_df.rename(
            columns={
                "Project Objective": "objective",
                "Cell Line of Choice": "cell_line",
                "Project Scope": "scope",
            },
            inplace=True,
        )
    except Exception as e:
        print(e)

    try:
        bill_df.rename(columns={"Scope": "scope",}, inplace=True)
    except Exception as e:
        print(e)

    bill_df = bill_df.reset_index(drop=True)
    # print(bill_df)
    # print(bill_df.columns.values)
    return bill_df


def gen_pi_df(pi_name, bill_df):
    pi_df = bill_df[bill_df["PI"] == pi_name]
    proj_list = pi_df["proj_no"].unique().tolist()
    new_proj_list = []
    for proj in proj_list:
        new_proj_list.append(proj.split(" ")[0])
    return pi_df, new_proj_list


def create_email_obj(email_body, email_signature, email_subject, email_recipients):
    outlook = win32.Dispatch("outlook.application")
    email_obj = outlook.CreateItem(0)
    # print()
    # print(email_recipients)
    # print(set(email_recipients))
    # print()
    email_obj.To = "; ".join(set(email_recipients))
    # mail.cc = recipient2
    email_obj.bcc = "Porter, Shaina"
    email_obj.Subject = email_subject
    email_obj.HtmlBody = "<br>".join([email_body, email_signature])
    email_obj.Display(False)

    return email_obj


def fetch_attachment(proj_no):
    try:
        path = CORE_PROJECTS_DIR
        for name in glob.glob(os.path.join(path, "*{}".format(proj_no))):
            folder = name
        # print("CORE Project folder is called {}".format(os.path.basename(folder)))
        os.chdir(folder)
        ppt_list = glob.glob("*.pptx")
        ppt_list = [ppt for ppt in ppt_list if "~$" not in ppt]
        latest_ppt = os.path.join(folder, max(ppt_list, key=os.path.getctime))

        # print("Slidedeck path = {}".format(latest_ppt))
        return latest_ppt
    except:
        print("couldn't find slidedeck in CORE Project folder")
        print("Project Number = {}".format(proj_no))
        return None


def add_attachments(mail_obj, proj_list):
    for proj_no in proj_list:
        attachment = fetch_attachment(proj_no)
        if attachment is not None:
            try:
                mail_obj.Attachments.Add(str(attachment))
            except Exception as e:
                print(f"Could not add attachment: {attachment}")
                print(e)

    return mail_obj


def add_scope_info():
    return


def email_single_project(pi_df, signature):
    print("  This is a single-project email.")

    def text_fxn(user, PI, gene, scope):
        PI_name = PI.split(", ")[1]
        user_name = user.split(", ")[1]
        salutations = f"Hi {PI_name}"
        if user.strip() != PI.strip():
            salutations += f" and {user_name}"
        if scope == "Edited cell pool":
            text = (
                f"{salutations},<br><br>Great news! Your "
                + gene
                + " edited cell pool project is complete and ready for pickup.  Please see the attached slide deck for details.<br><br>The last slide is the most informative.  We were able to get over XX% total editing in the pool with ~XX% out of frame indels.<br><br>We have a contactless pickup system in place right now.  Please coordinate with XXXXX to let XXXX know a good time window for you to pick up these cells. During the agreed upon time, XXX will place the frozen vials of cells in a dry ice bucket in M4170. The bucket will be on the counter in front of you when you walk in.  Just push the door to open.  It is always unlocked.  If you would like the live cultures as well, please come in the next day or so.  The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.<br><br>Don't hesitate to contact me if you have any questions.<br>Best,<br><br>SM"
            )
        elif scope == "Cell Fitness/Dependency Assay":
            text = f"{salutations},<br><br>Great news! Your {gene} fitness assay is complete. Please see the attached slide deck for details.<br><br>We do/do not see a strong dependency for this gene in this background.<br><br>Please let me know if you have any questions.<br><br>Best,<br>SM"
        elif scope == "gRNA with validation":
            text = f"{salutations},<br><br>Attached are the designs, off-target analysis, and our recommendations for which gRNAs to move forward with for your {gene} {scope} project.<br><br>Please let me know if you have any questions or if you would like to move forward with our recommendations.<br><br>Best,<br>SM<br>"
        else:
            text = (
                f"{salutations},<br><br>Great news! Your "
                + gene
                + " project is complete and ready for pick up.  Please see the attached slide deck for details.<br> We currently have a contactless pickup system in place.  Please arrange a time window with XXXXX in which someone can pick up the cells.  At the agreed upon time, he/she will place your frozen vials of cells into a dry ice bucket in M4170.  The dry ice bucket will be on the counter in front of you as you walk in.  Your live cultures will be in the first incubator to the right (top incubator, bottom shelf) and labeled accordingly. Please also bring dry ice for the pickup.<br><br>As always, please let me know if you have any questions.<br><br>Best,<br>SM"
            )
        return text

    def cell_line_text(user, PI, cell_line, gene, objective, scope):
        PI_name = PI.split(", ")[1]
        user_name = user.split(", ")[1]
        salutations = f"Hi {PI_name}"
        if user.strip() != PI.strip():
            salutations += f" and {user_name}"
        if scope == "Edited cell pool":
            text = f"""
            {salutations},<br>
            <br>
            Great news! Your {gene} {cell_line} edited cell pool project is complete and ready for pick up.  Please see the attached slide deck for details.<br>
            <br>
            The last slide is the most informative.  We were able to get over XX% total editing in the pool with ~XX% out of frame indels.<br>
            <br>
            We have a contactless pickup system in place right now.  Please coordinate with XXXXX to let XXXX know a good time window for you to pick up these cells. During the agreed upon time, XXXX will place the frozen vials of cells in a dry ice bucket in M4170. The bucket will be on the counter in front of you when you walk in.  The door is always unlocked.  If you would like the live cultures as well, please come in the next day or so.  The live cultures will be in the incubator to the right as you walk in (top incubator, bottom shelf).  Please bring dry ice for the pickup.<br>
            <br>
            Don't hesitate to contact me if you have any questions.<br>
            Best,<br>
            <br>
            SM
            """
        elif scope == "Cell Fitness/Dependency Assay":
            text = f"""{salutations},<br>
            <br>
            Great news! Your {gene} {cell_line} fitness assay is complete. Please see the attached slide deck for details.<br>
            <br>
            We do/do not see a strong dependency for this gene in this background.<br>
            <br>
            Please let me know if you have any questions.<br>
            <br>
            Best,<br>
            SM
            """
        else:
            text = f"""
            {salutations},<br>
            <br>
            Great news! Your {gene} {objective} project in {cell_line} cells is complete and ready for pick up.  Please see the attached slide deck for details.<br>
            We currently have a contactless pickup system in place.  Please arrange a time window with XXXXX in which someone can pick up the cells.  At the agreed upon time, he/she will place your frozen vials of cells into a dry ice bucket in M4170.  The dry ice bucket will be on the counter in front of you as you walk in.  Your live cultures will be in the first incubator to the right (top incubator, bottom shelf) and labeled accordingly. Please also bring dry ice for the pickup.<br><br>As always, please let me know if you have any questions.<br>
            <br>
            Best,<br>
            SM
            """
        return text

    def Emailer(text, signature, subject, recipient, recipient2, proj_no):
        outlook = win32.Dispatch("outlook.application")
        mail = outlook.CreateItem(0)
        if recipient == recipient2:
            mail.To = recipient
        else:
            mail.To = recipient + ";" + recipient2
        # mail.cc = recipient2
        mail.bcc = "Porter, Shaina"
        mail.Subject = subject
        mail.HtmlBody = "<br>".join([text, signature])
        mail.Display(False)
        # try:
        #     path = CORE_PROJECTS_DIR
        #     for name in glob.glob(os.path.join(path, "*{}".format(proj_no))):
        #         folder = name
        #     print("CORE Project folder is called {}".format(os.path.basename(folder)))
        #     os.chdir(folder)
        #     ppt_list = glob.glob("*.pptx")
        #     latest_ppt = os.path.join(folder, max(ppt_list, key=os.path.getctime))
        #     path1 = latest_ppt
        #     print("Slidedeck path = {}".format(path1))
        #     attachment = path1
        #     mail.Attachments.Add(attachment)
        # except:
        #     print("couldn't find slidedeck in CORE Project folder")
        #     print("Project Number = {}".format(proj_no))
        mail = add_attachments(mail, [proj_no])

    for i in range(pi_df["SRM_no"].count()):
        first_row = pi_df.head(1)
        # print(first_row)
        # print(first_row.columns)
        SRM_no = first_row["SRM_no"].values[0]
        PI = first_row["PI"].values[0]
        user = first_row["user"].values[0]
        gene = first_row["gene"].values[0]
        scope = first_row["scope"].values[0]
        try:
            proj_no = first_row["proj_no"].values[0].split(" ")[0]
        except:
            print(
                "Please check that your project numbers are entered correctly in the CAGEServices_Excel Export"
            )
            sys.exit()
        print(f"  * SRM# {SRM_no}  ({proj_no})")
        # print(cell_line+gene+objective)
        if scope in [
            "Cell Line Creation",
            "Edited cell pool",
            "Cell Fitness/Dependency Assay",
        ]:
            try:
                cell_line = first_row["cell_line"]
                objective = first_row["objective"]
                text = cell_line_text(user, PI, cell_line, gene, objective, scope)
            except:
                text = text_fxn(user, PI, gene, scope)
        else:
            text = text_fxn(user, PI, gene, scope)

        if scope == "Edited cell pool":
            subject = gene + " edited cell pool complete"
        elif scope == "Cell Line Creation":
            subject = gene + " cell line complete"
        elif scope == "Cell Fitness/Dependency Assay":
            subject = f"Cell Fitness/Dependency Assay for {gene}"
        elif scope == "gRNA with validation":
            subject = f"gRNA Designs for {gene}, SRM order {str(SRM_no)}"
        else:
            subject = f"{gene} Project Complete"
        Emailer(str(text), signature, subject, user, PI, proj_no)
    return


def email_multi_project(pi_df, signature, proj_list):
    print("  This is a multi-project email.")

    def group_by_scope():
        scope_dfs = []
        for scope in pi_df["scope"].unique().tolist():
            scope_dfs.append(pi_df[pi_df["scope"] == scope])
        return scope_dfs

    def bullet_point_text(row):
        bullet_point = f"<li>SRM {row['SRM_no']}: {row['gene']} ({row['scope']})</li>"
        print(f"  * SRM # {row['SRM_no']}  ({row['proj_no']})")
        return bullet_point

    def text_fxn(recipients_string, list_of_bullet_points):
        text = f"""
            Hi {recipients_string},<br>
            <br>
            Great news! The following projects are now complete and ready for pickup for your lab. Please see the attached slide decks for details.<br>
            <ul>
            {''.join(list_of_bullet_points)}
            </ul>
            We currently have a contactless pickup system in place.  Please arrange a time window with XXXXX in which someone can pick up the cells.  At the agreed upon time, he/she will place your frozen vials of cells into a dry ice bucket in M4170.  The dry ice bucket will be on your right as you walk in.  Your live cultures will be in the top incubator next to the dry ice bucket on the bottom shelf. Please also bring dry ice for the pickup.<br>
            <br>
            As always, please let me know if you have any questions.<br>
            <br>
            Best,<br>
            SM
            """
        return text

    def grna_text(recipients_string, list_of_bullet_points):
        text = f"""
            Hi {recipients_string},<br>
            <br>
            Attached are the designs, off-target analysis, and our recommendations for which gRNAs to move forward with for the following projects:<br>
            <ul>
            {''.join(list_of_bullet_points)}
            </ul>
            Please let me know if you have any questions or if you would like to move forward with our recommendations.<br>
            <br>
            Best,<br>
            SM
            """
        return text

    PI = pi_df.head(1)["PI"].values[0]
    pi_first = PI.split(", ")[1]
    recipients_string = pi_first

    user_names = pi_df["user"].unique().tolist()
    # print(user_names)
    if len(user_names) > 0:
        try:
            user_firsts = [name.split(", ")[1] for name in user_names]
        except:
            user_firsts = user_names

        for i, name in enumerate(user_firsts):
            if len(user_firsts) > 1:
                recipients_string += ", "
            if i == len(user_firsts) - 1:
                if len(user_firsts) == 1:
                    recipients_string += " "

                recipients_string += "and "
            recipients_string += name
    # print(user_firsts)

    # print(type(PI))
    # print(PI)
    # print(pi_df[["SRM_no"]])
    # print(pi_df[["gene"]])
    # print(pi_df[["scope"]])
    # print(pi_df.columns)
    list_of_bullet_points = pi_df.apply(bullet_point_text, axis=1)
    list_of_recipients = set([PI, *pi_df["user"].unique().tolist()])

    # all_scopes = pi_df["scope"].unique().tolist()
    # if len(all_scopes) == 1 and all_scopes[0] == "gRNA with validation":
    if "gRNA Approval Status" in pi_df.columns:
        email_subject = "gRNA Designs Ready"
        email_body = grna_text(recipients_string, list_of_bullet_points)
    else:
        if len(list_of_bullet_points) == 1:
            email_subject = "Project Ready for Pickup"
        else:
            email_subject = "Projects Ready for Pickup"
        email_body = text_fxn(recipients_string, list_of_bullet_points)

    email_obj = create_email_obj(
        email_body, signature, email_subject, list_of_recipients
    )
    email_obj = add_attachments(email_obj, proj_list)
    return email_obj


def handle_pi_df(pi_df, signature, proj_list):
    # if len(pi_df) > 1:
    #     email_multi_project(pi_df, signature, proj_list)
    #     return
    # else:
    #     email_single_project(pi_df, signature)
    #     return

    email_multi_project(pi_df, signature, proj_list)
    return


def main():
    signature = """<br>
        <br>
        Shondra Miller, Ph.D. <br>
        Director, Center for Advanced Genome Engineering (CAGE)<br>
        Associate Director, Shared Resources for Comprehensive Cancer Center<br>
        Associate Member, Department of Cell and Molecular Biology<br>
        St. Jude Children's Research Hospital<br>
        262 Danny Thomas Place <br>
        Memphis, TN  38105<br>
        Office: (901)595-7313 - Rm M4413<br>
        Lab:  (901)595-7314 – Rm M4160<br>
        """
    bill_df = gen_bill_df()
    print()
    for pi_name in bill_df["PI"].unique().tolist():
        print(f"Working on email for {pi_name}")
        pi_df, proj_list = gen_pi_df(pi_name, bill_df)

        handle_pi_df(pi_df, signature, proj_list)
        print()

    print("All draft emails are created")
    return


if __name__ == "__main__":
    main()
