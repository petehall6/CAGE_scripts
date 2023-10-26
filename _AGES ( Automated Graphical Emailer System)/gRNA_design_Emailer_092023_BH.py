"""

04-23-2020
@author: Rachel Levine
Program to create emails for sending out gRNA Design slidedeck to PI and requestor
Requires input of excel export from SRM2 from Approval page of CAGE PROJECTS

"""

from __future__ import print_function
import pandas as pd
import numpy as np
import os
import glob
from datetime import date
import win32com.client as win32

os.chdir(os.path.dirname(os.path.abspath(__file__)))
list_of_files = glob.glob("CAGEServices_Excel Export*.xls")

latest_file = max(list_of_files, key=os.path.getctime)
NGS_df = pd.read_excel(latest_file, sheet_name = 0) #opens first sheet of excel file
NGS_df.to_csv("email.csv", encoding='utf-8-sig')
NGS_df = pd.read_csv("email.csv")

NGS_df.rename(columns={
	"SRM Project #": "SRM_no",
	"Project Number": "proj_no",
	"Requested By": "user",
	"Scope": "scope",
	"Species": "species",
	"Target Gene Name": "gene",
	"gRNA Design - Final PowerPoint File": "ppt"
	}, inplace=True)

NGS_df = NGS_df.reset_index(drop=True)

def Emailer(text, subject, recipient,recipient2,proj_no,ppt):

	signature = "<br><br>Baranda Hansen, MS<br>Center for Advanced Genome Engineering (CAGE)<br>St. Jude Children's Research Hospital<br>Department of Cell and Molecular Biology<br>262 Danny Thomas Place<br>Advanced Research Center M4160<br>Memphis, TN  38105<br>Phone: (901)-595-7314"
	outlook = win32.Dispatch('outlook.application')
	mail = outlook.CreateItem(0)
	mail.To = recipient
	mail.cc = recipient2
	mail.Subject = subject
	mail.HtmlBody = text+signature
	mail.Display(False)
	try:
		path = r"Z:\\ResearchHome\\Groups\\millergrp\\home\\common\\CORE PROJECTS"

		for name in glob.glob(os.path.join(path,"*{}".format(proj_no))):
			folder = name
		print("CORE Project folder is called {}".format(os.path.basename(folder)))
		path1=folder+"/"+ppt
		print("Slidedeck Location = {}".format(path1))
		attachment  = path1
		mail.Attachments.Add(attachment)
	except:
		print ("couldn't find slidedeck in CORE Project folder")
    
def text_fxn(PI,user,gene,scope):#number, gene, edit, PI, results, proj_no):
	PI_name = PI.split(', ')[1]
	user_name = user.split(', ')[1]
	#signature = "<br><br>Rachel M Levine, Ph.D.<br>Center for Advanced Genome Engineering (CAGE)<br>St. Jude Children's Research Hospital<br>Department of Cell and Molecular Biology<br>262 Danny Thomas Place<br>Memphis, TN  38105<br>Office: (901)595-7314 - Rm D5052"
	text = "Hi "+ PI_name +" and "+user_name+",<br><br>Attached are the designs, off-target analysis, and our recommendations for which gRNAs to move forward with for your "+gene+" "+scope+" project.<br><br>Please let me know if you have any questions or if you would like to move forward with our recommendations.<br><br>Best,<br>SM<br>"
	return text

for i in range(NGS_df['SRM_no'].count()):
	SRM_no = int(NGS_df.loc[i,'SRM_no'])
	proj_no = NGS_df.loc[i,'proj_no'].split(" ")[0]
	PI = NGS_df.loc[i,'PI']
	user = NGS_df.loc[i,'user']
	scope = NGS_df.loc[i,'scope']
	gene = NGS_df.loc[i,'gene']
	ppt = NGS_df.loc[i,'ppt']
	text = text_fxn(PI,user,gene,scope)
	Emailer(str(text), 'gRNA Designs for '+gene+', SRM order '+str(SRM_no),PI,user,proj_no,ppt)

print("All draft emails are created")