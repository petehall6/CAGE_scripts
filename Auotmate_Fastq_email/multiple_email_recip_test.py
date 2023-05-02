import os
import pandas as pd
import shutil
import win32com.client

#way to scroll through NGS run date, CAGE#, well# and pull fastq files for automated zipping and emailing.

NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"
DEST_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/Python/CAGE_Programs/fastq_emailer/zip_output"

#get input
run_date = input("Enter run date: ")
CAGE_NUM = input("Enter project Number: ")
email_recip = input("Recipient email address: ")
email_cc = input("CC'd email address")

def _get_and_zip():
    working_dir = os.path.join(NGS_DIR,run_date,"joined")
    os.chdir(working_dir)

    #finds all subdirectories in working_dir and if cage_num in the name appends to cagenum_dirs list.  List comprehensions don't need 'and' for multiple conditions
    print("Searching for project folder.  This may take a moment.")
    #TODO add exception handling if no file is found
    cagenum_dirs = [ dir.name for dir in os.scandir(working_dir) if dir.is_dir() if CAGE_NUM in dir.name]
    
    #create a df from cagenum_dirs for easier reading
    temp_df = pd.DataFrame(cagenum_dirs, columns=['Directory'])
    #set first row to display as 1
    dirs_df = temp_df
    dirs_df.index +=1
    print()
    print(dirs_df)
    print()
    #index still initilazes to 0. Subtract 1 to account for print/actual format
    #TODO
    dir_choice = int(input("Enter directory choice number: "))-1
    print()
    cage_dir = dirs_df.iloc[dir_choice,0]

    print(f"CAGE Directory: {cage_dir}")

    #point back to NGS folder to scrape fastq file names from CSV
    target_dir = os.path.join(working_dir,cage_dir)
    os.chdir(target_dir)


    #finds all .csv's in target dir by looking at file extension. Doesn't need to be a list but it's easier to work with.
    csv_list = [ f.name for f in os.scandir(target_dir) if f.name.endswith(".csv") ]


    #reads CSV into a df to parse plate names
    csv_df = pd.read_csv(os.path.join(target_dir,csv_list[0]))
    #converts plate names into list
    fastq_names = csv_df['Name'].values.tolist()
    print()
    print(f"FASTQ target names: {fastq_names}")

    #points back to joined folder and loops through fastq_names list to see if fastq files name match.  Startswith requies str or tuple. Tuple will iterate through loop by itself.
    print("Finding FASTQ files")
    fastq_zips = [f for f in os.scandir(working_dir) if f.name.startswith(tuple(fastq_names))]

    #point to destination folder and zip files in fastq_zips list
    os.chdir(DEST_DIR)
    #creates folder of files to be zipped
    to_zip_dir = CAGE_NUM+"_fastq"
    zipped_name = CAGE_NUM+"_fastq_zipped"
    
    os.mkdir(to_zip_dir)
    
    print("Copying files over.")
    for fastq in fastq_zips:
        shutil.copy(fastq, to_zip_dir)

    print("Copy sucessfull")
    print("Zipping files")
    shutil.make_archive(zipped_name, "zip", to_zip_dir)
    print("Zipping complete")

def _emailer(email_recip):
    
    body = f"""
    Hi,
    
    I have attached the requested FASTQ files for project{CAGE_NUM}.  If you have any questions please feel free to contact me.  Thanks.
        
    """

    outlook = win32com.client.Dispatch("Outlook.Application")
    email = outlook.CreateItem(0)
    email.To = email_recip
    email.CC = email_cc
    email.Subject = "FASTQ Emailer Test"
    email.Body = body
    #explicitly converting slashes because I'm lazy right now
    attachment = os.path.join(DEST_DIR,CAGE_NUM+"_fastq_zipped.zip").replace("\\","/")
    print(attachment)
    email.Attachments.Add(attachment)
    email.Display(True)
    #email.Save()
    #email.Send()
    
def _delete_after():
    #shutil.rmtree since its a not an empty directory
    shutil.rmtree(os.path.join(DEST_DIR,CAGE_NUM+"_fastq").replace("\\","/"))
    #os.remove because rmtree doesnt work on .zip directory but os.remove does...go figure
    os.remove(os.path.join(DEST_DIR,CAGE_NUM+"_fastq_zipped.zip").replace("\\","/"))
    
_get_and_zip()

_emailer(email_recip)

_delete_after()

