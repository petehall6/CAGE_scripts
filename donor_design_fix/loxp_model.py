from __future__ import print_function
from Bio.SeqUtils import GC
from ssODN_for_deletions import add_features
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
from Bio.SeqUtils import nt_search
from operator import itemgetter, attrgetter, methodcaller
import glob
import copy
import datetime
import sys
import re




def display_IDT(name, full_gene_record, gRNA_cut_site, HA_length, len_loxP_RE):
    print(
        "Copy below to paste into IDT website, CHOOSE bulk input and CSV Format in top option"
    )
    ssODN_len = len(
        full_gene_record.seq[
            gRNA_cut_site - HA_length : gRNA_cut_site + len_loxP_RE + HA_length
        ]
    )
    sense = str(
        full_gene_record.seq[
            gRNA_cut_site - HA_length : gRNA_cut_site + len_loxP_RE + HA_length
        ]
    )
    antisense = str(
        full_gene_record.seq[
            gRNA_cut_site - HA_length : gRNA_cut_site + len_loxP_RE + HA_length
        ].reverse_complement()
    )
    # insert_locs=[ssODN_len-1, ssODN_len-2,ssODN_len+2-ssODN_len, ssODN_len+1-ssODN_len]
    # for star in insert_locs:
    #    sense=sense[:star]+'*'+sense[star:]
    #    antisense=antisense[:star]+'*'+antisense[star:]
    print(str("{}.sense.ssODN,{}".format(name, sense)))
    print(str("{}.antisense.ssODN,{}\n".format(name, antisense)))
    print(str("{}.sense.ssODN".format(name)))
    print(str("{}.antisense.ssODN\n".format(name)))




def add_primers_to_ssODN_record(ssODN_gene_record, current_gbk, name, folder_name1):
    # Go through the current gbk, pull out primer sequences "primer bind" and add to ssODN_gene_record
    # primer_dict={}
    primer_list = []
    for feature in current_gbk.features:
        if feature.type == "primer_bind":
            primer_list.append(feature)
            if feature.extract(str(current_gbk.seq)) in ssODN_gene_record:
                primer = feature.extract(str(current_gbk.seq))
                try:
                    ssODN_gene_record.features.append(
                        SeqFeature(
                            FeatureLocation(
                                ssODN_gene_record.seq.find(primer),
                                ssODN_gene_record.seq.find(primer) + len(primer),
                            ),
                            type="primer_bind",
                            strand=1,
                            qualifiers={
                                "label": feature.qualifiers["label"][0],
                                "note": feature.qualifiers["label"][0],
                            },
                        )
                    )
                except KeyError:
                    ssODN_gene_record.features.append(
                        SeqFeature(
                            FeatureLocation(
                                ssODN_gene_record.seq.find(primer),
                                ssODN_gene_record.seq.find(primer) + len(primer),
                            ),
                            type="primer_bind",
                            strand=1,
                            qualifiers={
                                "label": feature.qualifiers["note"][0],
                                "note": feature.qualifiers["note"][0],
                            },
                        )
                    )
            elif (
                feature.extract(current_gbk.seq).reverse_complement()
                in ssODN_gene_record
            ):
                #print("found the reverse complement:{}".format(feature))
                try:
                    primer = feature.extract(current_gbk.seq).reverse_complement()
                    ssODN_gene_record.features.append(
                        SeqFeature(
                            FeatureLocation(
                                ssODN_gene_record.seq.find(primer),
                                ssODN_gene_record.seq.find(primer) + len(primer),
                            ),
                            type="primer_bind",
                            strand=-1,
                            qualifiers={
                                "label": feature.qualifiers["label"][0],
                                "note": feature.qualifiers["label"][0],
                            },
                        )
                    )
                except KeyError:
                    ssODN_gene_record.features.append(
                        SeqFeature(
                            FeatureLocation(
                                ssODN_gene_record.seq.find(primer),
                                ssODN_gene_record.seq.find(primer) + len(primer),
                            ),
                            type="primer_bind",
                            strand=-1,
                            qualifiers={
                                "label": feature.qualifiers["note"][0],
                                "note": feature.qualifiers["note"][0],
                            },
                        )
                    )
            else:
                print("not found")
        else:
            pass
        ssODN_gene_record.annotations={"molecule_type": "DNA"}
        filepath = str("/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/"+folder_name1+'/')
        directory = str(os.path.join(filepath, f"{name}.loxP.sense.ssODN.gbk"))
        output_file = open(directory, 'w')
        SeqIO.write(ssODN_gene_record, directory,"genbank")
        output_file.close()

        
def RE_site(loxP, current_gbk, feature):
    restriction_enzyme_dict = {
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "EcoRI": " GAATTC",
        "NotI": "GCGGCCGC",
        "EcoRI": "GAATTC",
    }
    RE_search_region = current_gbk.seq[
        feature.location.start - 300 : feature.location.start + 330
    ]
    for (
        enzyme,
        site,
    ) in (
        restriction_enzyme_dict.items()
    ):  # Check 300bp upstream/downstream for a unique restriction enzyme
        if site not in RE_search_region:
            RE_pick = enzyme
            RE_seq = site
            RE_side = str(input("Is this loxP on 5' or 3'? (Enter 5 or 3)->")).strip()
            if RE_side == "5":
                insertion_dict = {
                    "len_1": len(loxP),
                    "label_1": "loxP",
                    "len_2": len(site),
                    "label_2": enzyme,
                }
                loxP_RE = loxP + site
                break
            elif RE_side == "3":
                insertion_dict = {
                    "len_1": len(site),
                    "label_1": enzyme,
                    "len_2": len(loxP),
                    "label_2": "loxP",
                }
                loxP_RE = site + loxP
            else:
                print("Not a valid option")
            break
        else:
            print(
                "Warning-no unique restriction site found +/-300 bp (enzymes searched: {}".format(
                    str(restriction_enzyme_dict.keys())
                )
            )
            pass
    return loxP_RE, insertion_dict






def add_HA_and_cut(full_gene_record,gRNA_cut_site,HA_length,insertion_dict,loxP_RE,name):
    
    full_gene_record.features.append(
        SeqFeature(
            FeatureLocation(gRNA_cut_site - HA_length, gRNA_cut_site),
            type="misc_feature",
            qualifiers={"label": "Left_HA", "note": "Left_HA"},
        )
    )
    full_gene_record.features.append(
        SeqFeature(
            FeatureLocation(gRNA_cut_site, gRNA_cut_site + insertion_dict["len_1"]),
            type="misc_feature",
            qualifiers={
                "label": insertion_dict["label_1"],
                "note": insertion_dict["label_1"],
            },
        )
    )
    full_gene_record.features.append(
        SeqFeature(
            FeatureLocation(
                gRNA_cut_site + insertion_dict["len_1"],
                gRNA_cut_site + insertion_dict["len_1"] + insertion_dict["len_2"],
            ),
            type="misc_feature",
            qualifiers={
                "label": insertion_dict["label_2"],
                "note": insertion_dict["label_2"],
            },
        )
    )
    full_gene_record.features.append(
        SeqFeature(
            FeatureLocation(
                gRNA_cut_site + len(loxP_RE), gRNA_cut_site + len(loxP_RE) + HA_length
            ),
            type="misc_feature",
            qualifiers={"label": "Right_HA", "note": "Right_HA"},
        )
    )
    full_gene_record.features.append(
        SeqFeature(
            FeatureLocation(
                gRNA_cut_site - HA_length, gRNA_cut_site + len(loxP_RE) + HA_length
            ),
            type="misc_feature",
            qualifiers={"label": name + ".sense.ssODN", "note": "Full_ssODN"},
        )
    )
    return full_gene_record





def get_user_input_gRNA(current_gbk):
    # Lists gRNAs in the project mod_NGS.gbk and queries user for choice
    gRNA_dict = {}  # Dict of found gRNAs:features
    for feature in current_gbk.features:
        if feature.type == "misc_feature":
            pat = "[A-Z]+\d+[.][A-Z0-9]+([.]g\d+)"  # Pattern to recongize Project_number.gene.gRNA. (i.e. CAGE52.hWASSUP.g1)
            pat2 = "[A-Z]+\d+\.[^\.]+\.([^\.]+)\.(g\d+)"   # Pattern to recongize Project_number.gene.gRNA_type.gRNA. (i.e. CAGE52.hWASSUP.9.g1  or CAGE52.hWASSUP.12a.g1)
            try:
                name = str(feature.qualifiers["label"]).strip("['']")
                gRNA_pat = re.match(pat, name, flags=re.I)
                if gRNA_pat ==None:
                    gRNA_pat = re.match(pat2, name, flags=re.I)
            except KeyError:
                name = str(feature.qualifiers["note"]).strip("['']")
                feature.qualifiers["label"] = [name]  # Add gRNA as label if it does not exist already to make upstream functions easier
                gRNA_pat = re.match(pat, name, flags=re.I)
                if gRNA_pat ==None:
                    gRNA_pat = re.match(pat2, name, flags=re.I)
            if gRNA_pat != None:
                gRNA_dict[name] = feature
    while True:
        # Menu for picking a gRNA.  Somewhat complicated b/c has number choice to pick gRNA rather than writing out full name
        temp_dict = {}
        for k, v in enumerate(gRNA_dict):
            print("{} : {}".format(k + 1, v))
            temp_dict[str(k + 1)] = v
        gRNA = input("Enter target gRNA for loxp site integration # (q to go back) -> ")
        if str(gRNA) in temp_dict.keys():
            return (
                temp_dict[str(gRNA)],
                gRNA_dict[temp_dict[str(gRNA)]],
            )  # Return gRNA name and gRNA feature location
        elif gRNA.upper() == "Q":
            return "Q", "NA"
            break
        else:
            print("Not a valid choice")


def create_loxP_from_core_projects_dir(filepath, folder_name1):
    
    current_gbk = SeqIO.read(filepath, "genbank")

    print("\n\nEnter first gRNA for loxP site integration\n")
    gRNA1, feature1 = get_user_input_gRNA(
    current_gbk
    )  
    gRNA_cut_site_dict = {
        -1: 6,
        1: 17,
    }  # Position of gRNA cut sites, -1 =nontemplate strand, 1=template strang
    
    gRNA_cut_site1 = feature1.location.start + gRNA_cut_site_dict[feature1.location.strand]

    print("\n\nEnter second gRNA for loxP site integration\n")
    gRNA2, feature2 = get_user_input_gRNA(
    current_gbk
    )
    gRNA_cut_site2 = feature2.location.start + gRNA_cut_site_dict[feature2.location.strand]

    loxP1 = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
    loxP2 = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"

    print(f"for the first gRNA: {gRNA1} ")
    loxP_RE1, insertion_dict1 = RE_site(loxP1, current_gbk, feature1)
    print(f"for the second gRNA: {gRNA2} ")
    loxP_RE2, insertion_dict2 = RE_site(loxP2, current_gbk, feature2)


    
    if abs(gRNA_cut_site1 - gRNA_cut_site2) <=40:
        print("gRNAs too close")
        return

    if gRNA_cut_site1 < gRNA_cut_site2:
        ssODN_gene_record = (
        current_gbk.seq[:gRNA_cut_site1] + loxP_RE1 + current_gbk.seq[gRNA_cut_site1:gRNA_cut_site2] + loxP_RE2 + current_gbk.seq[gRNA_cut_site2:]
    )
    else:
        ssODN_gene_record = (
        current_gbk.seq[:gRNA_cut_site2] + loxP_RE2 + current_gbk.seq[gRNA_cut_site2:gRNA_cut_site1] + loxP_RE1 + current_gbk.seq[gRNA_cut_site1:]
    )

    HA_length = 40 

    full_gene_record = SeqRecord(seq=ssODN_gene_record, id="current_gbk.id")
    full_gene_record.annotations["topology"] = "linear"


    name = str(gRNA1)
    add_HA_and_cut(full_gene_record,gRNA_cut_site1,HA_length,insertion_dict1,loxP_RE1,name)
    name = str(gRNA2)
    if gRNA_cut_site1 < gRNA_cut_site2:
        gRNA_cut_site2=gRNA_cut_site2 + len(loxP_RE1)
    
    
    add_HA_and_cut(full_gene_record,gRNA_cut_site2,HA_length,insertion_dict2,loxP_RE2,name)
    
    #print(full_gene_record)


    if gRNA1 != "Q" and gRNA2 != "Q":
        add_primers_to_ssODN_record(
            full_gene_record, current_gbk, name, folder_name1
        )

    print("\nSelected gRNAs below:\n")
    print(gRNA1, feature1)
    print(gRNA2, feature2)


    

    





def lox_p_Main():
    print("\n*Note* new genbank will be created from first project ID, gRNAs from the second project ID will be added\n")
    CAGEnum1 = input("Enter your *first project ID (case sensitive): ")
    userfile1, folder_name1, filename1 = findfile(CAGEnum1)

    CAGEnum2 = input("Enter your *second project ID (case sensitive): ")
    userfile2, folder_name2, filename2 = findfile(CAGEnum2)

    if userfile1 == -1:
        return
    if userfile2 == -1:
        return

    gRNAs1, gRNA_loc1, gRNA_strand1 = gRNA_Finder(userfile1)
    gRNAs2, gRNA_loc2, gRNA_strand2 = gRNA_Finder(userfile2)
    #print("\ngrnaloc2 ", gRNA_loc2)

    genome=SeqIO.read(userfile1,'genbank')
    genome2 = SeqIO.read(userfile2,'genbank')
    gRNA_filename = str(filename1)+'_+_'+str(filename2)+'.all_gRNAs'
    filepath = add_new_features(genome, gRNAs2, gRNA_loc2, gRNA_strand2, gRNA_filename, folder_name1, genome2)

    gRNA1, gRNA2, g1_loc, g2_loc, g1_strand, g2_strand = gRNA_Finder_for_deletions(filepath)
    gRNA1 = gRNA1[0]
    gRNA2 = gRNA2[0]
    

    ssODN_and_filename = str(gRNA1)+'_'+str(gRNA2)+'.del.sense.ssODN'
    print('\n\n\n Name: \n\n\n', ssODN_and_filename)
    genome=SeqIO.read(filepath,'genbank') 
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
    #print(g1_start, g1_end)
    #print(g2_start, g2_end)

    if g1_strand == 1:
        g1_cut =int(g1_start) + 17
    if g1_strand == -1:
        g1_cut = int(g1_start) + 6
    if g2_strand == 1:
        g2_cut = int(g2_start) + 17
    if g2_strand == -1:
        g2_cut = int(g2_start) + 6
    print('g1 strand is: {}, and cut will be: {}'.format(g1_strand,g1_cut))
    print('g2 strand is: {}, and cut will be: {}'.format(g2_strand,g2_cut))
    #print(g1_cut, g2_cut)

    if g1_cut < g2_cut:
        ssODN_slice1 = genome.seq[(g1_cut-40):g1_cut]
        ssODN_slice2 = genome.seq[g2_cut:(g2_cut+40)]

        genome_slice1 = genome.seq[(0):g1_cut]
        genome_slice2 = genome.seq[g2_cut:(len(genome)-1)]

        #print(int(len(genome.seq)-1))
        first=1
        margin = (g2_cut-g1_cut)

    if g1_cut > g2_cut:
        ssODN_slice1 = genome.seq[(g2_cut-40):g2_cut]
        ssODN_slice2 = genome.seq[g1_cut:(g1_cut+40)]


        genome_slice1 = genome.seq[(0):g1_cut]
        genome_slice2 = genome.seq[g2_cut:(len(genome)-1)]

        #print(int(len(genome.seq)-1))
        first=2
        margin = (g1_cut-g2_cut)

    #bridging_ssODN = genome_slice1 + '   JOINING SITE   '+ genome_slice2
    bridging_ssODN = ssODN_slice1 + ssODN_slice2
    new_genome = genome_slice1 + genome_slice2
    #create deletion genbank
    add_features(bridging_ssODN, genome, g1_cut, g2_cut, filepath, new_genome, filename1, ssODN_and_filename, ssODN_slice1, ssODN_slice2, folder_name1, first)
    
    create_loxP_from_core_projects_dir(filepath,folder_name1)



    #iterate through existing features, adding them to new deletion genbank as long as they are not in the deleted sequence
def add_new_features(genome, gRNAs2, gRNA_loc2, gRNA_strand2, ssODN_and_filename, folder_name1, genome2):
    #print('\n\n\n\n\ngrnas 2 list',gRNAs2, '\n\n\n\n\n')
    sequence_string = str(genome.seq[0:(len(genome)-1)])
    # Create a sequence
    
    sequence_object = Seq(sequence_string)
    # Create a record
    record = SeqRecord(
        sequence_object,
        id= genome.id,
        name=genome.name,
        annotations={"molecule_type": "DNA"},
        description=genome.description )#+ '\na deletion between' + gRNA1 + 'and'+ gRNA2+ 'in'+ filename1)

    
    filepath = str("/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/"+folder_name1+'/'+ssODN_and_filename+'_mark.gbk')
    output_file = open(filepath, 'w')
    SeqIO.write(record, filepath, 'genbank')
    output_file.close()
    new_genome = SeqIO.read(filepath,'genbank')
    for feat in genome.features:
        
        feat_loc = feat.location
        feat_loc = str(feat_loc)
        #print("\n\nfeature location: ", feat_loc, '\n\n', "feature: ", feat)
        if "," in feat_loc:
            continue
        start, end = feat_loc.split(":")
        start = start.strip("[")
        start = start.strip("<")
        start = start.strip(">")

        end = end.strip("] (+)")
        end = end.strip("] (-)")
        end = end.strip("<")
        end = end.strip(">")
        feat_seq = genome.seq[int(start):int(end)]

        align_count = genome.seq.count(feat_seq, start=0)
        start = new_genome.seq.find(feat_seq, start=0)
        end = start + len(feat_seq)
        reg_feat_index = genome.seq.find(feat_seq, start=0)
        if align_count > 1:
            #print("\nsequence was found more than once, Feature:\n ", feat)
            continue
                
        rev_feat = feat_seq.reverse_complement()
        align_count_rev = genome.seq.count(rev_feat, start=0)
        if align_count_rev > 1:
            print("\nreverse complement of sequence was found more than once, Feature:\n ", rev_feat, feat)
            continue
                
        new_feature = SeqFeature(FeatureLocation(int(start), int(end)), id = feat.id, strand = feat.strand, type = feat.type, qualifiers=feat.qualifiers)                 
        new_genome.features.append(new_feature)
        #print("adding updated feature:", new_feature)

    #add grnas from second genbank file
    for i in range(1,len(gRNAs2)+1):
        print(i)
    
        gRNA2 = gRNAs2.get(int(i))
        print("gRNA2",gRNA2)
        gRNA2_site = str(gRNA_loc2.get(int(i)))
        g2_strand = gRNA_strand2.get(int(i))
        print("grna2 site",gRNA2_site)

        start, end = gRNA2_site.split(":")
        start = start.strip("[")
        end = end.strip("] (+)")
        end = end.strip("] (-)")
        
        gRNA_seq = genome2.seq[int(start):int(end)]
        g_index = new_genome.seq.find(gRNA_seq, start=0)

        rev_comp = gRNA_seq.reverse_complement()
        g_rev_index = new_genome.seq.find(rev_comp, start=0)

        if g_index == -1 and g_rev_index == -1:
            print(gRNA2,"not found")
            continue

        if g_rev_index == -1:
            direction = "left"
            grna_feature = SeqFeature(FeatureLocation(int(g_index), int(g_index) + len(gRNA_seq)), id = gRNA2, type='misc_feature', qualifiers= {'label': gRNA2, 'Sequence': gRNA_seq, 'direction': direction})

        else:
            direction = "right"
            grna_feature = SeqFeature(FeatureLocation(int(g_index), int(g_index) + len(gRNA_seq)), id = gRNA2, strand = g2_strand, type='misc_feature', qualifiers= {'label': gRNA2, 'Sequence': gRNA_seq, 'direction': direction})

        new_genome.features.append(grna_feature)

    gbk_file = open(filepath, "w")
    SeqIO.write(new_genome, filepath, "genbank")
    gbk_file.close()
    print("genbank file created: ", filepath)


    return filepath




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
            print("\n\nFolder Found:",folder,",   looking through folder...")
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
                    print(' ',count, ":", filename)
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
    #print(userfile)
    #print(folder_name)
    #print(filename1)
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
            #print(qualifierDict)
            
            for key in qualifierDict:
                if key == 'label' or key == 'note':
                    #print(qualifierDict.get(key))
                    if ('.g') in str(qualifierDict.get(key)):
                        #print('\n', qualifierDict.get(key))
                        
                        Dict = {count:qualifierDict.get(key)}
                        gRNAs.update(Dict)
                        #print(count, ":", qualifierDict.get(key))

                        gRNA_loc = {count: feat.location}
                        gRNA_location.update(gRNA_loc)
                        #print("guide",qualifierDict.get(key), "\nlocation",gRNA_loc)
                        

                        gRNA_strand = {count: feat.strand}
                        gRNA_strand.update(gRNA_strand)

                        count = count + 1
    if count == 1:
        print("no gRNAs found!!")
        quit()


    return gRNAs, gRNA_location, gRNA_strand

  

def gRNA_Finder_for_deletions(filepath):
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
           
            for key in qualifierDict:
                if key == 'label' or key == 'note':
                    #print(qualifierDict.get(key))
                    if ('.g') in str(qualifierDict.get(key)):
                        #print('\n', qualifierDict.get(key))
                        
                        Dict = {count:qualifierDict.get(key)}
                        gRNAs.update(Dict)
                        print(count, ":", qualifierDict.get(key))

                        dict2 = {count: feat.location}
                        gRNA_location.update(dict2)

                        dict3 = {count: feat.strand}
                        gRNA_strand.update(dict3)

                        count = count + 1
    if count == 1:
        print("no gRNAs found!!")
        quit()
    else:
        userNum1 = input("Select your first gRNA (for deletion) (type number and return (ex. 2)): ")
        #if str(userNum1) not in gRNAs.keys:
            #print('number entered not associated with a gRNA')
            #return
        gRNA1 = gRNAs.get(int(userNum1))
        g1_loc = gRNA_location.get(int(userNum1))
        g1_strand = gRNA_strand.get(int(userNum1))
        userNum2 = input("Select your second gRNA (for deletion) (type number and return (ex. 2)): ")
        #if str(userNum2) not in gRNAs.keys:
            #print('number entered not associated with a gRNA')
            #return
        gRNA2 = gRNAs.get(int(userNum2))
        g2_loc = gRNA_location.get(int(userNum2))
        g2_strand = gRNA_strand.get(int(userNum2))


    return gRNA1, gRNA2, g1_loc, g2_loc, g1_strand, g2_strand

 

if __name__ == "__lox_p_Main__":
    lox_p_Main()