#pull all primers made for gene and add to gbk file

import pandas as pd
from Bio import SeqFeature, SeqRecord, SeqIO
from Bio import GenBank
from Bio.GenBank import Record
from Bio.SeqFeature import FeatureLocation
from Bio import Entrez
import os
from pathlib import PureWindowsPath
import shutil

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




target_gbk = gene+"_mod_NGS.gbk"
grna_seqs= []
target_files=[]
target_frames = []
target_projects = []
found_primers = []
found_guides = []
final_features = []
print(f"File Target: {target_gbk}")
print("Scanning for folders.")

#create target list project
for f in os.scandir(core_projects_dir):
    if f.is_dir() and gene in f.name:
        target_projects.append(f.name)

print(f"Number of associated projects: {len(target_projects)}")
print("List of projects: "+str(target_projects))


#open each gbk file and create dataframe.  parse grna and add to grna list
for dir in target_projects:
    os.chdir(os.path.join(core_projects_dir,dir))
    for f in os.scandir(os.path.join(core_projects_dir,dir)):
            if f.name == target_gbk:
                target_files.append(f)
            #TODO permissions not given at this moment
            #if f.name.endswith('.gbk'):
                #shutil.copyfile(f, dest)
print(f"Number of target gbks: {len(target_files)}")
print(str(target_files))



#open gbk file and find primer_bind feature:
gbk = SeqIO.read(target_files[1], "genbank")

#look for primer_bind feature

for feature in gbk.features:
    #find primers
    if feature.type =="primer_bind":
        found_primers.append(feature)
    #if feature.type == "misc_feature" and gene in feature.qualifiers

#TODOprint(f"Number of primers found: {len(found_primers)}")
#add features to genebank file
#add primers and guides separately
for feat in found_primers:
    
    #mark the location
    start_pos = feat.location.start
    end_pos = feat.location.end
    my_feature_location = FeatureLocation(start_pos, end_pos)
    
    #mark type
    my_feature_type = "primer_bind"
    
    #create SeqFeature
    my_feature = SeqFeature.SeqFeature(my_feature_location, type=my_feature_type)   
    
    final_features.append(my_feature)

#explicitly point back to custom folder for saving file
#os.chdir(dest)

#append final features to new gbk
#TODOprint(final_features)

#Load 'gene.gbk' from a core projects folder and print it
gbk_template_dir = os.path.join(core_projects_dir,dir)
gbk_template = gbk_template_dir+"\\gene.gbk"

template_gb_record = SeqIO.read(open(gbk_template, "r"), "genbank")
    
template_gb_record.features.append(final_features)
template_gb_record.annotations["molecule_type"] = "DNA"   


print(template_gb_record)

with open('all_primers.gbk', 'w') as output_file:
    SeqIO.write(template_gb_record, output_file, "genbank")





print("Completed")
