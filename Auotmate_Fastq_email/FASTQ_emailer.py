#5/4/2023 PMH

#This script compiles FASTQ files, zips them and attaches them to an email




import os
import pandas as pd
import shutil
import win32com.client

#way to scroll through NGS run date, CAGE#, well# and pull fastq files for automated zipping and emailing.
NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"
DEST_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/Python/CAGE_Programs/fastq_emailer/zip_output"

def _input():
    #get input
    run_date = input("Enter run date: ").strip()
    working_dir = os.path.join(NGS_DIR, run_date, "joined").replace("\\","/")

    #check if run_date is correct/folder exsists
    while os.path.isdir(working_dir) == True:
        print("Run date found.")
        break
        
    else:
        print("Run date not found.  Please reenter run date.")
        run_date = input("Enter run date: ").strip()
        working_dir = os.path.join(NGS_DIR, run_date, "joined").replace("\\","/")
        
    
    cagenum_dirs =[]
    #check project number/cage number.  If there anything in the cagenum/projectnum dirs list loop breaks and contiues to emailer
    while len(cagenum_dirs) == 0:
        cage_num = input("Please enter project number: ").upper().strip()
        print("Searching for project folder.  This may take a moment.")
        cagenum_dirs = [dir.name for dir in os.scandir(working_dir) if dir.is_dir() if cage_num in dir.name]
        print(f"Number of projects found: {len(cagenum_dirs)}")
        #empty list evaluate to false
        if not cagenum_dirs:
            print(f"Did not find project number {cage_num}.")    
            
    return working_dir,cagenum_dirs,cage_num

def _get_fastq_and_zip(working_dir, cagenum_dirs, cage_num):
    #create a df from cagenum_dirs for easier reading
    temp_df = pd.DataFrame(cagenum_dirs, columns=['Directory'])
    #set first row to display as 1
    dirs_df = temp_df
    dirs_df.index +=1
    print()
    print(dirs_df)
    print()
    #index still initilazes to 0. Subtract 1 to account for print/actual format
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
    to_zip_dir = cage_num+"_fastq"
    zipped_name = cage_num+"_fastq_zipped"
    
    #deletes dir if already exsits
    if os.path.exists(to_zip_dir):
        shutil.rmtree(to_zip_dir)
    os.mkdir(to_zip_dir)
    
    print("Copying files over.")
    for fastq in fastq_zips:
        shutil.copy(fastq, to_zip_dir)

    #using .make_archive to zip the folder and only folder/contents.
    #just using zipfile will result in zip bomb (entire file path is zipped)
    print("Copy sucessfull")
    print("Zipping files")
    shutil.make_archive(zipped_name, "zip", to_zip_dir)
    print("Zipping complete")
    
#opens outlook and emails primary(To:) aka email_recip and secondary (CC:)
def _emailer(cage_num):
    email_recip = input("Recipient email address: ").strip().replace(" ","").split(",")
    email_cc = input("CC'd email address: ").strip().replace(" ","").split(",")
    
    body = (
    f"""
    Hi,
    
    I have attached the requested FASTQ files for project {cage_num}.  If you have any questions please feel free to contact me.  Thanks.    
    """
    
    )

    outlook = win32com.client.Dispatch("Outlook.Application")
    email = outlook.CreateItem(0)
    #splits each element from list and seperates each with ';' allowing for multiple recipenants and cc's
    email.To = ";".join(email_recip)
    email.CC = ";".join(email_cc)
    email.Subject = f"{cage_num} FASTQ Files"
    email.Body = body
    #explicitly converting slashes for ease
    attachment = os.path.join(DEST_DIR,cage_num+"_fastq_zipped.zip").replace("\\","/")
    email.Attachments.Add(attachment)
    email.Display(True)
    #uncomment line below to automatically send email without opening outlook window
    #email.Send()
    
#remove old zipped files because cant be trusted to clean up after themselves
def _delete_after(cage_num):
    #shutil.rmtree since its a not an empty directory
    shutil.rmtree(os.path.join(DEST_DIR,cage_num+"_fastq").replace("\\","/"))
    #os.remove because rmtree doesnt work on .zip directory but os.remove does...go figure
    os.remove(os.path.join(DEST_DIR,cage_num+"_fastq_zipped.zip").replace("\\","/"))


working_dir, cagenum_dirs,cage_num = _input()
        
_get_fastq_and_zip(working_dir, cagenum_dirs, cage_num)

_emailer(cage_num)


clean = input("Do you want to delete zipped files? (y/n): ").upper()

if clean != "N":
    _delete_after(cage_num)

print("Program completed.")

