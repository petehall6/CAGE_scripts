#pull all primers made for gene and add to gbk file


from Bio import SeqIO
import os
from pathlib import PureWindowsPath
import shutil
import glob

core_projects_dir = PureWindowsPath(r"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS")
dest = os.path.dirname(os.path.abspath(__file__))
os.chdir(core_projects_dir)


""" if species == "H":
    gene = "h"+raw_gene
elif species =="M":
    gene = "m"+raw_gene

print(f"Target gene: {gene}") """


#template_gbk = os.path.join(dest,f"{gene}_copied.gbk")
target_gbk = "_mod_NGS.gbk"
target_files=[]
target_projects = []
found_features = []
copied_template = False
#print(f"File Target: {target_gbk}\n\n")
#print("Scanning for folders.\n\n")

#output_file = os.path.join(dest,"all_features.gbk")
#gene_search_string = f"{core_projects_dir}\\{gene}-*\\"
print("Finding Targets")


#create target list project by searching core project folders for gene name
for f in glob.glob("*"):
    if os.path.isdir(f):
        target_projects.append(f)





print(f"Number of associated project folders: {len(target_projects)}")


#create list of the target
for dir in target_projects:
    os.chdir(os.path.join(core_projects_dir,dir))
    #print(os.path.join(core_projects_dir, dir))
    for f in os.scandir(os.path.join(core_projects_dir,dir)):
            if f.name.endswith('NGS.gbk'):
                target_files.append(f)

print(f"Target files found: {len(target_files)}")
#point back to the source folder to pull template gbk

os.chdir(dest)


for file in target_files:
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
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

print(f"\nNumber of features found: {len(found_features)}\n")



print("Completed\n")

