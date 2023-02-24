#pull all primers made for gene and add to gbk file


from Bio import SeqIO
import os
from pathlib import PureWindowsPath
import shutil
import glob

core_projects_dir = PureWindowsPath(r"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS")
dest = os.path.dirname(os.path.abspath(__file__))
os.chdir(core_projects_dir)


raw_gene = input("Enter target gene: ").upper()
species = input("Enter (h)uman or (m)ouse: ").upper()

if species == "H":
    gene = "h"+raw_gene
elif species =="M":
    gene = "m"+raw_gene

print(f"Target gene: {gene}")


template_gbk = os.path.join(dest,f"{gene}_copied.gbk")
target_gbk = raw_gene+"_mod_NGS.gbk"
target_files=[]
target_projects = []
found_features = []
copied_template = False
print(f"File Target: {target_gbk}\n\n")
print("Scanning for folders.\n\n")

output_file = os.path.join(dest,"all_features.gbk")
gene_search_string = f"{core_projects_dir}\\{gene}-*\\"



#create target list project by searching core project folders for gene name
target_projects= glob.glob(gene_search_string)

#found at least one project with associated gene
if len(target_projects) > 1:
    print(f"Number of associated projects: {len(target_projects)}")
    print(f"List of projects: {target_projects}\n")

#create list of the target
    for dir in target_projects:
        os.chdir(os.path.join(core_projects_dir,dir))
        for f in os.scandir(os.path.join(core_projects_dir,dir)):
                if f.name == target_gbk:
                    target_files.append(f)
                if f.name == "gene.gbk" and f.name.endswith('.gbk') and copied_template is False: #only copies the gbk from the first project for speed
                    shutil.copy(f, os.path.join(dest,f"{gene}_copied.gbk"))
                    copied_template = True
    print("Template gbk created\n")

#point back to the source folder to pull template gbk
    os.chdir(dest)

    output_handle = open(f"{gene}_all_features.gbk", "w")

    fail_num = 0
#
    for file in target_files:
        for record in SeqIO.parse(file, "genbank"):
            for feature in record.features:
                if feature.type == "primer_bind":
                    found_features.append(feature)
            
            for feature in record.features:
                if feature.type == "CDS":
                    found_features.append(feature)

                if feature.type == "misc_feature":
                    try:
                        for qual in feature.qualifiers["note"]:
                            if "CAGE" in qual:
                                found_features.append(feature)
                    except KeyError as e:
                        failed_directory = os.path.abspath(file)
                        print(f"Failed to add feature from project @ : {failed_directory})")
                        failed_feature = str(feature.qualifiers)
                        print(f"No note added to feature: {failed_feature}")
                        fail_num = fail_num+1               
                        
    print(f"\nNumber of features not added: {fail_num}")
    print(f"\nNumber of features found: {len(found_features)}\n")


    for record in SeqIO.parse(template_gbk, "genbank"):
        record.features = found_features
        
    SeqIO.write(record, output_handle, "genbank")   

    output_handle.close()
    print("Removing copied files")
    os.remove(template_gbk)


    print("Completed\n")

else:
    print("Gene not found.  Rerun script.  Try to be more specific or double check spelling")