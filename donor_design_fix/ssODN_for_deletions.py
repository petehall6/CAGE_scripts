from Bio.SeqUtils import GC
from Bio import File, SeqIO
import os
import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from pandas.core.indexes.base import ensure_index_from_sequences
from off_target_all_RDM_v4 import off_target_all_primers
import shutil
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC


def Main():
    projectID = input("Enter your project ID (case sensitive): ")
    userfile, folder_name, filename1 = findfile(projectID)
    if userfile == -1:
        return
    gRNA1, gRNA2, g1_loc, g2_loc, g1_strand, g2_strand = gRNA_Finder(userfile)
    
    gRNA1 = gRNA1[0]
    gRNA2 = gRNA2[0]
    

    ssODN_and_filename = str(gRNA1)+'_'+str(gRNA2)+'.del.sense.ssODN'
    print('\n\n\n Name: \n\n\n', ssODN_and_filename)
    genome=SeqIO.read(userfile,'genbank') 
    #print(genome.features)
    print('gRNA1 strand: ', g1_strand, '\ngNRA2 strand: ', g2_strand)

    g1_loc = str(g1_loc)
    g2_loc = str(g2_loc)
    print("your selections are: ",gRNA1, gRNA2)
    print('location of', gRNA1,': ',g1_loc)
    print('location of', gRNA2,': ',g2_loc)
    

    g1_start, g1_end = g1_loc.split(":")
    g1_start = g1_start.strip("[")
    g1_end = g1_end.strip("] (+)")
    g1_end = g1_end.strip("] (-)")
    
    g2_start, g2_end = g2_loc.split(":")
    g2_start = g2_start.strip("[")
    g2_end = g2_end.strip("] (+)")
    g2_end = g2_end.strip("] (-)")
    print(g1_start, g1_end)
    print(g2_start, g2_end)

    
    #print("g1_strand",g1_strand)
    #print("g2_strand",g2_strand)
    if g1_strand == 1:
        g1_cut = int(g1_start) + 17
    if g1_strand == -1:
        g1_cut = int(g1_start) + 6
    if g2_strand == 1:
        g2_cut = int(g2_start) + 17
    if g2_strand == -1:
        g2_cut = int(g2_start) + 6

    print(g1_cut, g2_cut)

    if g1_cut < g2_cut:
        ssODN_slice1 = genome.seq[(g1_cut-40):g1_cut]
        ssODN_slice2 = genome.seq[g2_cut:(g2_cut+40)]

        genome_slice1 = genome.seq[(0):g1_cut]
        genome_slice2 = genome.seq[g2_cut:(len(genome)-1)]

        print(int(len(genome.seq)-1))
        first = 1

        margin = (g2_cut-g1_cut)

    
    if g1_cut > g2_cut:
        ssODN_slice1 = genome.seq[(g2_cut-40):g2_cut]
        ssODN_slice2 = genome.seq[g1_cut:(g1_cut+40)]

        print("ssODN slice 1:  ",ssODN_slice1)
        print("ssODN slice 2:  ",ssODN_slice2)
        genome_slice1 = genome.seq[(0):g2_cut]
        genome_slice2 = genome.seq[g1_cut:(len(genome)-1)]

        print(int(len(genome.seq)-1))
        first = 2

        margin = (g1_cut-g2_cut)

    
    
    #bridging_ssODN = genome_slice1 + '   JOINING SITE   '+ genome_slice2
    bridging_ssODN = ssODN_slice1 + ssODN_slice2
    
    new_genome = genome_slice1 + genome_slice2

    add_features(bridging_ssODN, genome, g1_cut, g2_cut, userfile, new_genome, filename1, ssODN_and_filename, ssODN_slice1, ssODN_slice2, folder_name, first)
    
    print("\nBridging ssODN: ", bridging_ssODN)
    #if first == 1:
        #print("cut site:  ",new_genome[g1_cut-40:g2_cut+40])
    #if first == 2:
        #print("cut site:  ",new_genome[g2_cut-40:g1_cut+40])
    

    
    
    #filepath = str("/rgs01/project_space/millergrp/CrisprDesigns/common/python/homology_arm_primer_finder/"+filename1+'_del.gbk')
    #del_genome=SeqIO.read(filepath,'genbank')

    #iterate through existing features, adding them to new deletion genbank as long as they are not in the deleted sequence
def add_features(bridging_ssODN, genome, g1_cut, g2_cut, original_userfile, new_genome_seq, filename, ssODN_name, ssODN_slice1, ssODN_slice2, folder_name, first):
    
    # Create a sequence
    sequence_string = str(new_genome_seq)
    sequence_object = Seq(sequence_string)
    # Create a record
    record = SeqRecord(
        sequence_object,
        id= genome.id,
        name=genome.name,
        annotations={"molecule_type": "DNA"},
        description=genome.description )#+ '\na deletion between' + gRNA1 + 'and'+ gRNA2+ 'in'+ filename1)

    #feature = SeqFeature(FeatureLocation(start=del_feat_index, end=int(del_feat_index + len(feat_seq))), type=feat.type, qualifiers=feat.qualifiers)
    #record.features.append(feature)
    
    # Save as GenBank file
    
    print("\n\n ssODN_name", ssODN_name)
    filepath = str("/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/"+folder_name+'/'+ssODN_name+'_mark.gbk')
    output_file = open(filepath, 'w')
    SeqIO.write(record, filepath, 'genbank')
    output_file.close()
    
    #print("\n\n 1 \n\n")
    new_genome = SeqIO.read(filepath,'genbank') 

    for feat in genome.features:
        #print("\n\n 2 \n\n")
        #print("original feature", feat)
        
        feat_loc = feat.location
        feat_loc = str(feat_loc)
        print("\n\nfeature location: ", feat_loc, '\n', "feature: ", feat)
        if "," in feat_loc:
            continue
        start, end = feat_loc.split(":")
        start = start.strip("[")
        end = end.strip("] (+)")
        end = end.strip("] (-)")
        feat_seq = genome.seq[int(start):int(end)]

        align_count = new_genome.seq.count(feat_seq, start=0)
        start = new_genome.seq.find(feat_seq, start=0)
        if start == -1:
            continue

        end = start + len(feat_seq)

        reg_feat_index = genome.seq.find(feat_seq, start=0)
        if align_count > 1:
            print("\nsequence was found more than once, Feature:\n ", feat)
            continue
        

        if first == 1:
            if int(start) >= g1_cut and int(start) <= g2_cut:
                print("can't add feature because it falls within deletion, Feature:\n", feat)
                continue
                
            if int(end) >= g1_cut and int(end) <= g2_cut:
                print("can't add feature because it falls within deletion, Feature:\n", feat)
                continue
        if first == 2:
            if int(start) >= g2_cut and int(start) <= g1_cut:
                print("can't add feature because it falls within deletion, Feature:\n", feat)
                continue
                
            if int(end) >= g2_cut and int(end) <= g1_cut:
                print("can't add feature because it falls within deletion, Feature:\n", feat)
                continue

            
                
        rev_feat = feat_seq.reverse_complement()
        align_count_rev = genome.seq.count(rev_feat, start=0)
        if align_count_rev > 1:
            print("\nreverse complement of sequence was found more than once, Feature:\n ", rev_feat, feat)
            continue
                
        new_feature = SeqFeature(FeatureLocation(int(start), int(end)), id = feat.id, strand = feat.strand, type = feat.type, qualifiers=feat.qualifiers)                 
        new_genome.features.append(new_feature)
        print("adding updated feature:", new_feature)

    
    ssODN_index = new_genome.seq.find(bridging_ssODN, start=0)
    ssODN_slice1_index = new_genome.seq.find(ssODN_slice1, start=0)
    ssODN_slice2_index = new_genome.seq.find(ssODN_slice2, start=0)

    if ssODN_index == -1:
        rev_comp = bridging_ssODN.reverse_complement()
        ssODN_index = new_genome.seq.find(rev_comp, start=0)

        ssODN_feature = SeqFeature(FeatureLocation(int(ssODN_index), int(ssODN_index) + len(bridging_ssODN)), id = ssODN_name, strand = -1, type='misc_feature', qualifiers= {'label': ssODN_name, 'Sequence': bridging_ssODN, 'direction': 'left'})
        ssODN_40bp_homology_arm_1 = SeqFeature(FeatureLocation(int(ssODN_slice1_index), int(ssODN_slice1_index) + len(ssODN_slice1)), id = '40 bp homology arm', strand = -1, type='misc_feature', qualifiers= {'label': '40 bp homology arm', 'Sequence': ssODN_slice1, 'direction': 'left'})
        ssODN_40bp_homology_arm_2 = SeqFeature(FeatureLocation(int(ssODN_slice2_index), int(ssODN_slice2_index) + len(ssODN_slice2)), id = '40 bp homology arm', strand = -1, type='misc_feature', qualifiers= {'label': '40 bp homology arm', 'Sequence': ssODN_slice2, 'direction': 'left'})
    ssODN_feature = SeqFeature(FeatureLocation(int(ssODN_index), int(ssODN_index) + len(bridging_ssODN)), id = ssODN_name, strand = 1, type='misc_feature', qualifiers= {'label': ssODN_name, 'Sequence': bridging_ssODN, 'direction': 'right'})
    ssODN_40bp_homology_arm_1 = SeqFeature(FeatureLocation(int(ssODN_slice1_index), int(ssODN_slice1_index) + len(ssODN_slice1)), id = '40 bp homology arm', strand = 1, type='misc_feature', qualifiers= {'label': '40 bp homology arm', 'Sequence': ssODN_slice1, 'direction': 'right'})
    ssODN_40bp_homology_arm_2 = SeqFeature(FeatureLocation(int(ssODN_slice2_index), int(ssODN_slice2_index) + len(ssODN_slice2)), id = '40 bp homology arm', strand = 1, type='misc_feature', qualifiers= {'label': '40 bp homology arm', 'Sequence': ssODN_slice2, 'direction': 'right'})
    
    
    
    new_genome.features.append(ssODN_feature)
    new_genome.features.append(ssODN_40bp_homology_arm_1)
    new_genome.features.append(ssODN_40bp_homology_arm_2)


    


    gbk_file = open(filepath, "w")
    SeqIO.write(new_genome, filepath, "genbank")
    gbk_file.close()
    print("genbank file created: ", filepath)

    return






#uses project ID to find correct folder in CORE PROJECTS, then enumerates possible file options in folder for the user to select from
#returns path to desired file and the name of the folder
def findfile(projectID):
    enumeratedFiles = {}
    enumeratedFolders = {}
    enumerated_filenames = {}
    count = 1
    #directory = "/research/groups/millergrp/home/common/CORE PROJECTS" #"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS"
    directory = "/research/groups/millergrp/home/common/CORE PROJECTS" #"/research/rgs01/project_space/millergrp/CrisprDesigns/common/python/homology_arm_primer_finder/CORE PROJECTS"

    #print(directory)
    #iterate through folders
    for folder in os.listdir(directory):
        #print(folder)
        #find project ID in folder name
        if projectID in folder: #if folder == projectID: #if projectID in folder: 
            #print("\n\nFolder Found:",folder,",   looking through folder...")
            #print(folder)
            #iterate through files in folder
            for filename in os.listdir(os.path.join(directory, folder)):
                filename_upper = filename.upper()
                #print(filename)
                #check for vector
                if filename_upper.endswith(".GBK"):
                    #print("\nfound",filename,"in",folder,"folder")
                    directory2 = os.path.join(directory, folder)
                    path = os.path.join(directory2, filename)
                    
                    Dict = {count: str(path)}
                    enumeratedFiles.update(Dict)

                    dict2 = {count: str(folder)}
                    enumeratedFolders.update(dict2)

                    dict3 = {count: str(filename)}
                    enumerated_filenames.update(dict3)

                    #print(enumeratedFiles)
                    print(count, ":", filename)
                    count = count + 1
                else:
                    continue
                  
    if count == 1:
        print("no files found!!")
        return -1, folder, -1
    else:
        userNum = input("\nWhich of the files above are you using? (type number and return (ex. 2)): ")
    #print(userNum)
    #print(enumeratedFiles)
    #print(enumeratedFiles.keys())
    #print(enumeratedFiles.get(int(userNum)))
    userfile = enumeratedFiles.get(int(userNum))
    folder_name = enumeratedFolders.get(int(userNum))
    filename1 = enumerated_filenames.get(int(userNum))
    print(userfile)
    print(folder_name)
    print(filename1)
    #print("Path to file:     ", userfile)
    return userfile, folder_name, filename1


# Uses filepath to open genbank file and find gRNA
def gRNA_Finder(filepath):
    genome=SeqIO.read(filepath,'genbank') 
    genomeLength = len(genome) - 1
    is_gRNA = False

    gRNAs = {}
    gRNA_location = {}
    gRNA_strand = {}
    count = 1
    
    for feat in genome.features:
        if feat.type == 'misc_feature':
            qualifierDict = feat.qualifiers
            #if qualifierDict.get('label') == None:
             #   continue
            for key in qualifierDict:
                if key == 'label' or key== "note":
                    #print(qualifierDict.get(key))
                    if ('.g') in str(qualifierDict.get(key)):
                        #print('\n', qualifierDict.get(key))
                        
                        Dict = {count:qualifierDict.get(key)}
                        gRNAs.update(Dict)
                        print(count, ":", qualifierDict.get(key)[0])

                        dict2 = {count: feat.location}
                        gRNA_location.update(dict2)

                        dict3 = {count: feat.strand}
                        gRNA_strand.update(dict3)

                        count = count + 1
    if count == 1:
        print("no gRNAs found!!")
        quit()
    else:
        userNum1 = input("Select your first gRNA (type number and return (ex. 2)): ")
        #if str(userNum1) not in gRNAs.keys:
            #print('number entered not associated with a gRNA')
            #return
        gRNA1 = gRNAs.get(int(userNum1))
        g1_loc = gRNA_location.get(int(userNum1))
        g1_strand = gRNA_strand.get(int(userNum1))
        userNum2 = input("Select your second gRNA (type number and return (ex. 2)): ")
        #if str(userNum2) not in gRNAs.keys:
            #print('number entered not associated with a gRNA')
            #return
        gRNA2 = gRNAs.get(int(userNum2))
        g2_loc = gRNA_location.get(int(userNum2))
        g2_strand = gRNA_strand.get(int(userNum2))


    return gRNA1, gRNA2, g1_loc, g2_loc, g1_strand, g2_strand

                        
if __name__ == "__Main__":
    Main()
            
