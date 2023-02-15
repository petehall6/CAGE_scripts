#pull all primers made for gene and add to gbk file

import pandas as pd
from Bio import SeqFeature, SeqRecord, SeqIO
from Bio import GenBank
from Bio.GenBank import Record
from Bio.SeqFeature import SeqFeature
import shutil

import os
from pathlib import PureWindowsPath

core_projects_dir = PureWindowsPath(r"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS")
dest = os.path.dirname(os.path.abspath(__file__))
os.chdir(core_projects_dir)

#TODOprint(os.getcwd())

for file in os.listdir(dest):
    if file.endswith(".gbk"):
        handle = file
#TODOprint(f"GBK handle: {handle}")
gene = "ADA2"
#gene = input("Enter target gene: *TEST* ADA2: ").upper()



template_gbk = "gene.gbk"
target_gbk = gene+"_mod_NGS.gbk"
target_files=[]
target_projects = []
found_features = []
print(f"File Target: {target_gbk}\n\n")
print("Scanning for folders.\n\n")

output_file = os.path.join(dest,"all_features.gbk")

#create target list project
for f in os.scandir(core_projects_dir):
    if f.is_dir() and gene in f.name:
        target_projects.append(f.name)

print(f"Number of associated projects: {len(target_projects)}")
print(f"List of projects: {target_projects}\n")


#open each gbk file and create dataframe.  parse grna and add to grna list
for dir in target_projects:
    os.chdir(os.path.join(core_projects_dir,dir))
    for f in os.scandir(os.path.join(core_projects_dir,dir)):
            if f.name == target_gbk:
                target_files.append(f)
            #TODO permissions not given at this moment
            #if f.name == "gene" and f.name.endswith('.gbk'):
                #shutil.copyfile(f, dest)
            
os.chdir(dest)

output_handle = open("all_features.gbk", "w")

for file in target_files:
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
            if feature.type == "primer_bind":
                found_features.append(feature)
        
        for feature in record.features:
            if feature.type == "CDS":
                found_features.append(feature)

            if feature.type == "misc_feature":
                for qual in feature.qualifiers["note"]:
                    if "CAGE" in qual:
                        found_features.append(feature)

    print(f"Features found: {len(found_features)}")


for record in SeqIO.parse(template_gbk, "genbank"):
    record.features = found_features
    
SeqIO.write(record, output_handle, "genbank")   




















""" for file in target_files:
    for record in SeqIO.parse(file, "genbank"):
        for feature in record.features:
            if feature.type == "primer_bind":
                found_features.append(feature)
                print(feature)
            
            if feature.type == "misc_feature":
                for qual in feature.qualifiers["note"]:
                    if "CAGE" in qual:
                        found_features.append(feature)
                        print(feature)

    print(f"Features found: {len(found_features)}")
    
    for record in SeqIO.parse(template_gbk, "genbank"):
        for feature in found_features:
            record.features.append(feature)
        
        SeqIO.write(record, output_handle, "genbank") """

output_handle.close()


print("Completed")
