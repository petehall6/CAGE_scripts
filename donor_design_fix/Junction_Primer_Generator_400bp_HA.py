from re import T
from Bio.SeqUtils import GC
from Bio import File, SeqIO
import os
import pandas as pd
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
from off_target_all_RDM_v4 import off_target_all_primers
import shutil
import numpy as np
# from json import dumps

SCRIPT_DIR = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs"


def jxn_primer_gen():
    # print("\n\n\n\nSTART:\n\n\n\n")
    names = []
    name_index_dict = {}
    is_reverse_complement = {}
    LHAprimers_Dict = {}
    RHAprimers_Dict = {}
    good_primers = {}
    projectID = input("Enter your project ID (case sensitive): ") 
    print()
    filepath, folder_name, filename = findfile(projectID)
    if filepath == -1:
        return
    source = filepath
    destination = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder_name+'/'+filename+'_backup_mark.gbk.BAK'
    
    if os.path.isfile(destination):
        print("  (Note: file backup already made.)")
    else:
        shutil.copyfile(source, destination)
    list = str(folder_name).split("-")
    ID = list[2]
    print()
    primer_special_name_YN = input("  Would you like to append an additional identifier to primer names? (Y or N): ")
    print()
    if primer_special_name_YN.upper() == 'Y':
        primer_ID = input("Please enter what you would like added to the primer names: ")
        ID = ID+'.'+primer_ID
        print()
    

    if filepath == -1:
        return
    else:
        LHAprimer3file, RHAprimer3file, leftHALen, rightHALen = homologyArmSeq(filepath, ID, folder_name)
        if LHAprimer3file == -1:
            return
        else:
            changeProductSize(leftHALen)
            LHAprimers_0, LHAprimers_1, LHAprimers_2, LHAprimers_3, LHAprimers_4, LHAprimers_5, LHAprimers_6, LHAprimers_7, LHAprimers_8, LHAprimers_9 = LHAprimer3call(ID, folder_name)

            #print(LHAprimers_Dict)
            changeProductSize(rightHALen)
            RHAprimers_0, RHAprimers_1, RHAprimers_2, RHAprimers_3, RHAprimers_4, RHAprimers_5, RHAprimers_6, RHAprimers_7, RHAprimers_8, RHAprimers_9 = RHAprimer3call(ID, folder_name)
            #print(RHAprimers_Dict)
            P5gen_F_min = 999999999999999999999999
            P5junc_DS_R_min = 999999999999999999999999
            P3junc_DS_F_min = 999999999999999999999999
            P3gen_R_min = 999999999999999999999999

            P5gen_F_best = ['1', '2', '3'] 
            P5junc_DS_R_best = ['1', '2', '3'] 
            P3junc_DS_F_best = ['1', '2', '3'] 
            P3gen_R_best = ['1', '2', '3']

            primer_sets = [   
                (LHAprimers_0, RHAprimers_0), 
                (LHAprimers_1, RHAprimers_1), 
                (LHAprimers_2, RHAprimers_2), 
                (LHAprimers_3, RHAprimers_3), 
                (LHAprimers_4, RHAprimers_4),
                (LHAprimers_5, RHAprimers_5), 
                (LHAprimers_6, RHAprimers_6), 
                (LHAprimers_7, RHAprimers_7), 
                (LHAprimers_8, RHAprimers_8), 
                (LHAprimers_9, RHAprimers_9)]

            CAGE_5_gen_F = ID+'.5gen.F'
            CAGE_5junc_DS_R = ID+'.5junc.DS.R'
            CAGE_3junc_DS_F = ID+'.3junc.DS.F'
            CAGE_3gen_R = ID+'.3gen.R'
            
            valid_primer_sets = []
            for (LHA_primers, RHA_primers) in primer_sets:
                lha_f_none = LHA_primers[CAGE_5_gen_F] == 'None'
                lha_r_none = LHA_primers[CAGE_5junc_DS_R] == 'None'
                rha_f_none = RHA_primers[CAGE_3junc_DS_F] == 'None'
                rha_r_none = RHA_primers[CAGE_3gen_R] == 'None'
            
                if any([lha_f_none, lha_r_none, rha_f_none, rha_r_none]):
                    pass
                else:
                    valid_primer_sets.append((LHA_primers, RHA_primers))
            primer_sets = valid_primer_sets
            
            print('primer sets are: {}'.format(primer_sets))
            primer_scores = [[999] * len(primer_sets), [999] * len(primer_sets)]
            
            for ind, (lha, rha) in enumerate(primer_sets):
                summary_frame = off_target_check_for_narrowing(filepath, lha, rha, folder_name)

                print(f"\n\nPrimer set {ind}:\n", summary_frame)
                print()

                P5gen_F_matches_2 = summary_frame.at[CAGE_5_gen_F,"2_mm"]
                P5junc_DS_R_matches_2 = summary_frame.at[CAGE_5junc_DS_R,"2_mm"]
                P3junc_DS_F_matches_2 = summary_frame.at[CAGE_3junc_DS_F,"2_mm"]
                P3gen_R_matches_2 = summary_frame.at[CAGE_3gen_R,"2_mm"]
                
                primer_scores[0][ind] = (
                    np.nansum([P5gen_F_matches_2, P5junc_DS_R_matches_2])
                ) / 2
                
                primer_scores[1][ind] = (
                    np.nansum([P3junc_DS_F_matches_2, P3gen_R_matches_2])
                ) / 2
                
            min_set5 = primer_scores[0].index(np.nanmin(primer_scores[0]))
            min_set3 = primer_scores[1].index(np.nanmin(primer_scores[1]))
            
            for lha, rha in [primer_sets[min_set5]]:

                summary_frame = off_target_check_for_narrowing(filepath, lha, rha, folder_name)

                P5gen_F_matches_0 = summary_frame.at[CAGE_5_gen_F,"0_mm"]
                P5junc_DS_R_matches_0 = summary_frame.at[CAGE_5junc_DS_R,"0_mm"]
                P3junc_DS_F_matches_0 = summary_frame.at[CAGE_3junc_DS_F,"0_mm"]
                P3gen_R_matches_0 = summary_frame.at[CAGE_3gen_R,"0_mm"]

                P5gen_F_matches_1 = summary_frame.at[CAGE_5_gen_F,"1_mm"]
                P5junc_DS_R_matches_1 = summary_frame.at[CAGE_5junc_DS_R,"1_mm"]
                P3junc_DS_F_matches_1 = summary_frame.at[CAGE_3junc_DS_F,"1_mm"]
                P3gen_R_matches_1 = summary_frame.at[CAGE_3gen_R,"1_mm"]

                P5gen_F_matches_2 = summary_frame.at[CAGE_5_gen_F,"2_mm"]
                P5junc_DS_R_matches_2 = summary_frame.at[CAGE_5junc_DS_R,"2_mm"]
                P3junc_DS_F_matches_2 = summary_frame.at[CAGE_3junc_DS_F,"2_mm"]
                P3gen_R_matches_2 = summary_frame.at[CAGE_3gen_R,"2_mm"]

                if (np.isnan(P5gen_F_matches_2)):
                    P5gen_F_min = 0
                    P5gen_F_best[0] = CAGE_5_gen_F
                    P5gen_F_best[1] = P5gen_F_matches_2
                    P5gen_F_best[2] = summary_frame.at[CAGE_5_gen_F,"sequence"]
                else:
                    P5gen_F_min = int(P5gen_F_matches_2)
                    P5gen_F_best[0] = CAGE_5_gen_F
                    P5gen_F_best[1] = P5gen_F_matches_2
                    P5gen_F_best[2] = summary_frame.at[CAGE_5_gen_F,"sequence"]


                if (np.isnan(P5junc_DS_R_matches_2)):
                    P5junc_DS_R_min = 0
                    P5junc_DS_R_best[0] = CAGE_5junc_DS_R
                    P5junc_DS_R_best[1] = P5junc_DS_R_matches_2
                    P5junc_DS_R_best[2] = summary_frame.at[CAGE_5junc_DS_R,"sequence"]
                else:
                    P5junc_DS_R_min = int(P5junc_DS_R_matches_2)
                    P5junc_DS_R_best[0] = CAGE_5junc_DS_R
                    P5junc_DS_R_best[1] = P5junc_DS_R_matches_2
                    P5junc_DS_R_best[2] = summary_frame.at[CAGE_5junc_DS_R,"sequence"]

            for lha, rha in [primer_sets[min_set3]]:

                summary_frame = off_target_check_for_narrowing(filepath, lha, rha, folder_name)

                P5gen_F_matches_0 = summary_frame.at[CAGE_5_gen_F,"0_mm"]
                P5junc_DS_R_matches_0 = summary_frame.at[CAGE_5junc_DS_R,"0_mm"]
                P3junc_DS_F_matches_0 = summary_frame.at[CAGE_3junc_DS_F,"0_mm"]
                P3gen_R_matches_0 = summary_frame.at[CAGE_3gen_R,"0_mm"]

                P5gen_F_matches_1 = summary_frame.at[CAGE_5_gen_F,"1_mm"]
                P5junc_DS_R_matches_1 = summary_frame.at[CAGE_5junc_DS_R,"1_mm"]
                P3junc_DS_F_matches_1 = summary_frame.at[CAGE_3junc_DS_F,"1_mm"]
                P3gen_R_matches_1 = summary_frame.at[CAGE_3gen_R,"1_mm"]

                P5gen_F_matches_2 = summary_frame.at[CAGE_5_gen_F,"2_mm"]
                P5junc_DS_R_matches_2 = summary_frame.at[CAGE_5junc_DS_R,"2_mm"]
                P3junc_DS_F_matches_2 = summary_frame.at[CAGE_3junc_DS_F,"2_mm"]
                P3gen_R_matches_2 = summary_frame.at[CAGE_3gen_R,"2_mm"]

                
                if (np.isnan(P3junc_DS_F_matches_2)):
                    P3junc_DS_F_min = 0
                    P3junc_DS_F_best[0] = CAGE_3junc_DS_F
                    P3junc_DS_F_best[1] = P3junc_DS_F_matches_2
                    P3junc_DS_F_best[2] = summary_frame.at[CAGE_3junc_DS_F,"sequence"]
                else:
                    P3junc_DS_F_min = int(P3junc_DS_F_matches_2)
                    P3junc_DS_F_best[0] = CAGE_3junc_DS_F
                    P3junc_DS_F_best[1] = P3junc_DS_F_matches_2
                    P3junc_DS_F_best[2] = summary_frame.at[CAGE_3junc_DS_F,"sequence"]

               
                if (np.isnan(P3gen_R_matches_2)):
                    P3gen_R_min = 0
                    P3gen_R_best[0] = CAGE_3gen_R
                    P3gen_R_best[1] = P3gen_R_matches_2
                    P3gen_R_best[2] = summary_frame.at[CAGE_3gen_R,"sequence"]
                else:
                    P3gen_R_min = int(P3gen_R_matches_2)
                    P3gen_R_best[0] = CAGE_3gen_R
                    P3gen_R_best[1] = P3gen_R_matches_2
                    P3gen_R_best[2] = summary_frame.at[CAGE_3gen_R,"sequence"]
                    

            print()
            print(f'Set with lowest avg. long_2:')
            print(f"\t5': Set {min_set5}")
            print(f"\t3': Set {min_set3}")
            
            print()
            
            #print(P5gen_F_best, P5junc_DS_R_best, P3junc_DS_F_best, P3gen_R_best)
            final_primers_df = pd.DataFrame([P5gen_F_best, P5junc_DS_R_best, P3junc_DS_F_best, P3gen_R_best], columns = ["Primer Name", "Number of matches for 2mm", "Primer Sequence"])
            print('\n\nFinal primers df:\n', final_primers_df.to_string(), '\n')
            to_continue = input("Above are the best primers next to the number of matches found for 2mm, Do you want this program to continue with these primers? (Y or N)")
            print()
            if to_continue.upper() == 'N':
                print('program aborting')
                return
            if to_continue.upper() != 'Y':
                print('Please type Y or N')


            LHAprimers_Dict[ID+".5gen.F"] = P5gen_F_best[2]
            LHAprimers_Dict[ID+".5junc.DS.R"] = P5junc_DS_R_best[2]

            RHAprimers_Dict[ID+".3junc.DS.F"] = P3junc_DS_F_best[2]
            RHAprimers_Dict[ID+".3gen.R"] = P3gen_R_best[2]

           
            genome=SeqIO.read(filepath,'genbank')
            
            
            for site, name in [
                (LHAprimers_Dict.get(ID+".5gen.F"), ID+".5gen.F"),
                (LHAprimers_Dict.get(ID+".5junc.DS.R"), ID+".5junc.DS.R"),
                (RHAprimers_Dict.get(ID+".3junc.DS.F"), ID+".3junc.DS.F"),
                (RHAprimers_Dict.get(ID+".3gen.R"), ID+".3gen.R")]:                
                index = 0
                align_count = genome.seq.count(site, start=index)
                if align_count > 1:
                    print("\n\n\n*** ALERT!! Primer sequence was found more than once in the genbank file ***\n\nname:",name,"\nseq:",site,"\nnumber of times found:",align_count)
                    return
                rev_primer = (Seq(site)).reverse_complement()
                align_count_rev = genome.seq.count(rev_primer, start=index)
                if align_count_rev > 1:
                    print("\n\n\n*** ALERT!! reverse complement of primer sequence was found more than once in the genbank file ***\n\nname:",name,"\nseq:",rev_primer,"\nnumber of times found:",align_count_rev)
                    return

                index = genome.seq.find(site, start=index)
                

                if index == -1:
                    print(f"{name} NOT FOUND! ...Trying reverse complement...")
                    Seq_primer = (Seq(site)).reverse_complement()
                    #print(Seq_primer)
                    index = 0
                    index = genome.seq.find(Seq_primer, start=index)
                    #print('\n\n', name, index, '\n\n')
                    
                    if index == -1:
                        print("*** Primer NOT FOUND!!! ***")
                        return
                    else:
                        print(f"{name} found at: {index}")
                        dict2 = {name: index}
                        name_index_dict.update(dict2)
                        rev_comp_dict = {name: True}
                        is_reverse_complement.update(rev_comp_dict)
                else:
                    print(f"{name} found at: {index}")
                    dict1 = {name: index}
                    name_index_dict.update(dict1)
            # Check if this program has added primers already, If true prints primers, runs GC content and off target analysis then aborts
            with open(filepath, "r") as gbk_file:
                    lines = gbk_file.readlines()
                    for line in lines:
                        if str('Sequence="'+RHAprimers_Dict.get(ID+".3gen.R")) in line:
                            #print(line)
                            print("\nThis program has already added primers to this gbk file!\n\n\n")
                            for feat in genome.features:
                            
                                qualifierDict = feat.qualifiers
                                if qualifierDict.get('label') == [ID+".5gen.F"] or qualifierDict.get('label') == [ID+".5gen.F2"] or qualifierDict.get('label') == [ID+".5gen.F3"]:
                                    name_5genF = qualifierDict.get('label')
                                    seq_5genF = qualifierDict.get('Sequence')
                                if qualifierDict.get('label') == [ID+".5junc.DS.R"] or qualifierDict.get('label') == [ID+".5junc.DS.R2"] or qualifierDict.get('label') == [ID+".5junc.DS.R3"]:
                                    name_5juncDSR = qualifierDict.get('label')
                                    seq_5juncDSR = qualifierDict.get('Sequence')
                                if qualifierDict.get('label') == [ID+".3junc.DS.F"] or qualifierDict.get('label') == [ID+".3junc.DS.F2"] or qualifierDict.get('label') == [ID+".3junc.DS.F3"]:
                                    name_3juncDSF = qualifierDict.get('label')
                                    seq_3juncDSF = qualifierDict.get('Sequence')
                                if qualifierDict.get('label') == [ID+".3gen.R"] or qualifierDict.get('label') == [ID+".3gen.R"] or qualifierDict.get('label') == [ID+".3gen.R"]:
                                    name_3genR = qualifierDict.get('label')
                                    seq_3genR = qualifierDict.get('Sequence')

                            gbk_file = open(filepath, "r")
                            #print(gbk_file.read())
                            gbk_file.close()
                            GC_Content(filepath, RHAprimers_Dict, LHAprimers_Dict, ID)
                            ismatch = off_target_check(filepath, LHAprimers_Dict,RHAprimers_Dict,folder_name)
                            if ismatch == -1:
                                return


                    
                            return
            #primers_df.to_clipboard(excel=True,index=False, header=False)

            # makes a list of all label names found in genbank file
            qualifier_label_list = []
            for feat in genome.features:
                if feat.type == 'primer_bind':
                    qualifierDict = feat.qualifiers
                    label = qualifierDict.get('label')
                    if label != None:
                        qualifier_label_list.append(label)

            #print('\n\n\nQUALIFIER LABEL LIST', qualifier_label_list)

            # If qualifier label list is empty, add name to names and create feature
            if not qualifier_label_list:
                print("List is empty")
                
                for site, name in [
                (LHAprimers_Dict.get(ID+".5gen.F"), ID+".5gen.F"),
                (LHAprimers_Dict.get(ID+".5junc.DS.R"), ID+".5junc.DS.R"),
                (RHAprimers_Dict.get(ID+".3junc.DS.F"), ID+".3junc.DS.F"),
                (RHAprimers_Dict.get(ID+".3gen.R"), ID+".3gen.R")]: 
                
                    names.append(name) 
                    final_index = str(name_index_dict.get(name))



                    if ("SSDNA" in filename.upper() or "MEGAMER" in filename.upper()) and name == ID+".5gen.F":
                        name = ID+".5gen.F"
                    if ("SSDNA" in filename.upper() or "MEGAMER" in filename.upper()) and name == ID+".3gen.R":
                        name = ID+".3gen.R"
                    #if (filename.upper().__contains__("SSDNA") or filename.upper.__contains__("MEGAMER")) and name == ID+".5gen.F":
                    #    name = ID+".5gen.F"
                    #if (filename.upper().__contains__("SSDNA") or filename.upper.__contains__("MEGAMER")) and name == ID+".3gen.R":
                    #    name = ID+".3gen.R"
                    if is_reverse_complement.get(name) == True:
                        new_feature = SeqFeature(FeatureLocation(int(final_index), int(final_index) + len(site)), type = 'primer_bind', id = name, strand = -1, qualifiers={'label': name, 'Sequence': site, 'direction': 'left'})
                    else:
                        new_feature = SeqFeature(FeatureLocation(int(final_index), int(final_index) + len(site)), type = 'primer_bind', id = name, strand = 1, qualifiers={'label': name, 'Sequence': site, 'direction': 'right'})
                    genome.features.append(new_feature)
                    print(new_feature)
                   

            # Looks for existing primers in file that were not added by this program
            # If found, adds incremented number to the end of primer name then creates feature and adds it to the genbank file
            else:
                for site, name in [
                (LHAprimers_Dict.get(ID+".5gen.F"), ID+".5gen.F"),
                (LHAprimers_Dict.get(ID+".5junc.DS.R"), ID+".5junc.DS.R"),
                (RHAprimers_Dict.get(ID+".3junc.DS.F"), ID+".3junc.DS.F"),
                (RHAprimers_Dict.get(ID+".3gen.R"), ID+".3gen.R")]:     
                    count = 0

                    for label in qualifier_label_list:
                        #print(label)
                        count = count + 1
                        if name in label:
                            lastchar = label[-1:]
                            if str(lastchar).isdigit():
                                lastchar = int(lastchar) + 1
                                new_name = name+str(lastchar)
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break     
                            elif name == label:
                                new_name = name+'2'
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break      

                        if str(ID+'.Junc.DS.R') in label and name == ID+'.5junc.DS.R':
                            lastchar = label[-1:]
                            if str(lastchar).isdigit():
                                lastchar = int(lastchar) + 1
                                new_name = name+str(lastchar)
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break
                            else:
                                new_name = name+'2'
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break      
                        if str(ID+'.Junc.DS.F') in label and name == ID+".3junc.DS.F":
                            lastchar = label[-1:]
                            if str(lastchar).isdigit():
                                lastchar = int(lastchar) + 1
                                new_name = name+str(lastchar)
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break
                            else:
                                new_name = name+'2'
                                name_index_dict[new_name] = name_index_dict.pop(name)
                                name = new_name
                                names.append(name)
                                break   
                    
                        #print('count:', count)
                        #print('length of qualifier label list:', len(qualifier_label_list))
                        if count == len(qualifier_label_list) and name not in names:
                            names.append(name)
                            break
                
                   
                    
                    final_index = str(name_index_dict.get(name))
                    ssDNA_list=['MEGAMER','SSDA']
                    if any(ext in filename.upper() for ext in ssDNA_list) and name == ID+".5gen.F":
                        name = ID+".5gen.F"
                    if any(ext in filename.upper() for ext in ssDNA_list) and name == ID+".3gen.R":
                        name = ID+".3gen.R"
                    #    )
                    #if (filename.upper().__contains__("SSDNA") or filename.upper.__contains__("MEGAMER")) and name == ID+".5gen.F":
                    #    name = ID+".5gen.F"
                    #if (filename.upper().__contains__("SSDNA") or filename.upper.__contains__("MEGAMER")) and name == ID+".3gen.R":
                    #    name = ID+".3gen.R"
                    #new_feature = SeqFeature(FeatureLocation(int(final_index), int(final_index) + len(site)), type = 'primer_bind', id = name, qualifiers={'label': name, 'Sequence': site})
                    if is_reverse_complement.get(name) == True:
                        new_feature = SeqFeature(FeatureLocation(int(final_index), int(final_index) + len(site)), type = 'primer_bind', id = name, strand = -1, qualifiers={'label': name, 'Sequence': site, 'direction': 'left'})
                    else:
                        new_feature = SeqFeature(FeatureLocation(int(final_index), int(final_index) + len(site)), type = 'primer_bind', id = name, strand = 1, qualifiers={'label': name, 'Sequence': site, 'direction': 'right'})
                    genome.features.append(new_feature)
                    print(new_feature)
            
            gbk_file = open(filepath, "w")
            SeqIO.write(genome, filepath, "genbank")
            gbk_file.close()
            gbk_file = open(filepath, "r")
            #print(gbk_file.read())
            gbk_file.close()
            GC_Content(filepath, RHAprimers_Dict, LHAprimers_Dict, ID)
            ismatch = off_target_check(filepath, LHAprimers_Dict,RHAprimers_Dict,folder_name)
            if ismatch == -1:
                return
            primers_table = [[names[0], LHAprimers_Dict.get(ID+".5gen.F")],
                                [names[1], LHAprimers_Dict.get(ID+".5junc.DS.R")],
                                [names[2], RHAprimers_Dict.get(ID+".3junc.DS.F")],
                                [names[3], RHAprimers_Dict.get(ID+".3gen.R")]]
            primers_df = pd.DataFrame(data=primers_table, columns=["Name", "Sequence"])
            #print(primers_df)
            #print(summary_frame)
            print('\nPrimers Table\n')

            if names[0] == ID+".5gen.F" and names[3] == ID+".3gen.R":
                print("\n ",LHAprimers_Dict.get(ID+".5gen.F"),': ', names[0], "\n ", LHAprimers_Dict.get(ID+".5junc.DS.R"),': ', names[1],"\n ",)
                print(" ",RHAprimers_Dict.get(ID+".3junc.DS.F"),': ', names[2],"\n ",RHAprimers_Dict.get(ID+".3gen.R"),': ', names[3], "\n ")
            
                print("\n ",LHAprimers_Dict.get(ID+".5gen.F"),"\n ",LHAprimers_Dict.get(ID+".5junc.DS.R"))
                print(" ",RHAprimers_Dict.get(ID+".3junc.DS.F"),"\n ",RHAprimers_Dict.get(ID+".3gen.R"),"\n")
                print(" ",names[0],"\n ",names[1])
                print(" ",names[2],"\n ",names[3],"\n\n")
                print("adding new primers\n\n\n")

            else:
                print("\n ",LHAprimers_Dict.get(ID+".5gen.F"),': ', names[0], "\n ", LHAprimers_Dict.get(ID+".5junc.DS.R"),': ', names[1],"\n ",)
                print(" ",RHAprimers_Dict.get(ID+".3junc.DS.F"),': ', names[2],"\n ",RHAprimers_Dict.get(ID+".3gen.R"),': ', names[3], "\n ")
            
                print("\n ",LHAprimers_Dict.get(ID+".5gen.F"),"\n ",LHAprimers_Dict.get(ID+".5junc.DS.R"))
                print(" ",RHAprimers_Dict.get(ID+".3junc.DS.F"),"\n ",RHAprimers_Dict.get(ID+".3gen.R"),"\n")
                print(" ",names[0],"\n ",names[1])
                print(" ",names[2],"\n ",names[3],"\n\n")
                print("adding new primers\n\n\n")
                

            return


# Runs a bowtie call for off target analysis
# writes a primer.fa file and moves it to the folder in core projects
# also moves bowtie files to the folder in core projects
def off_target_check(filepath, LHAprimers_Dict, RHAprimers_Dict, folder):
    print('primers are: {}, and {}'.format(LHAprimers_Dict, RHAprimers_Dict))
    LHAprimers_Dict.update(RHAprimers_Dict)
    #print("all primers:", LHAprimers_Dict)

    with open("primers.fa", "w") as f:
        for key, value in LHAprimers_Dict.items():
            f.write(str(">"+ key))
            f.write("\n")
            f.write(value)
            f.write("\n")
    for gb_record in SeqIO.parse(filepath,'genbank'):
        organism = gb_record.annotations['organism'].upper()
    ismatch = False
    print("organism:", organism)
    if 'SAPIENS' in organism:
        OTA = 'human'
        ismatch = True
    if 'MUSCULUS' in organism: 
        OTA = 'mouse'
        ismatch = True
    if 'RATTUS' in organism:
        OTA = 'rat'
        ismatch = True
    if 'SABAEUS' in organism:
        OTA = 'monkey'
        ismatch = True
    if ismatch == False:
        print('organism could not be found on genbank file, could not do off target analysis')
        quit()

    off_target_all_primers(os.getcwd(), OTA)

    return


def off_target_check_for_narrowing(filepath, LHAprimers_Dict, RHAprimers_Dict, folder):
    print('Primer dicdtionaries are: {} and {}'.format(LHAprimers_Dict, RHAprimers_Dict))
    LHAprimers_Dict.update(RHAprimers_Dict)
    print('after call')
    #print("all primers:", LHAprimers_Dict)
    print(LHAprimers_Dict)
    if 'None' in LHAprimers_Dict.values():
        print('NONE!!!')
    else:

        with open("primers.fa", "w") as f:
            for key, value in LHAprimers_Dict.items():
                f.write(str(">"+ key))
                f.write("\n")
                f.write(value)
                f.write("\n")
        for gb_record in SeqIO.parse(filepath,'genbank'):
            organism = gb_record.annotations['organism'].upper()
        ismatch = False
        # print("\norganism:", organism)
        if 'SAPIENS' in organism:
            OTA = 'human'
            ismatch = True
        if 'MUSCULUS' in organism: 
            OTA = 'mouse'
            ismatch = True
        if 'RATTUS' in organism:
            OTA = 'rat'
            ismatch = True
        if 'SABAEUS' in organism:
            OTA = 'monkey'
            ismatch = True
        if ismatch == False:
            print('organism could not be found on genbank file, could not do off target analysis')
            quit()

        #print(OTA)
        #print(organism)
        #f = open("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/primers.fa", "r")
        #print("primers.fa file here:", f.read)
        #f.close()
        summary_frame = off_target_all_primers(os.getcwd(), OTA)
        #source = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/primers.fa"
        #destination = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/primers_Mark.fa'
        #shutil.copyfile(source, destination)
        source1 = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/primer_0mm.bowtie"
        source2 ="/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/primer_1mm.bowtie"
        source3 = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/primer_2mm.bowtie"
        destination1 = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/primer_0mm_Mark.bowtie'
        destination2 = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/primer_1mm_Mark.bowtie'
        destination3 = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/primer_2mm_Mark.bowtie'
        #shutil.copyfile(source1, destination1)
        #shutil.copyfile(source2, destination2)
        #shutil.copyfile(source3, destination3)
        return summary_frame

# Finds GC content of the amplicons that would be created
def GC_Content(filepath, RHAprimers_Dict, LHAprimers_Dict, ID):
    genome=SeqIO.read(filepath,'genbank')
    leftHAleftPloc = ''
    
    for feat in genome.features: 
        qualifierDict = feat.qualifiers
        #print(LHAprimers_Dict)
        #print(ID+".5gen.F")
        #print(LHAprimers_Dict.get(ID+".5gen.F"))
        #print(LHAprimers_Dict.get(ID+".5junc.DS.R"))
        #print(RHAprimers_Dict.get(ID+".3junc.DS.F"))
        #print(RHAprimers_Dict.get(ID+".3gen.R"))
        # print("qualifier seq:", qualifierDict.get('Sequence'))
        
        if qualifierDict.get('Sequence') == [LHAprimers_Dict.get(ID+".5gen.F")]:
            leftHAleftPloc = str(feat.location)
        if qualifierDict.get('Sequence') == [LHAprimers_Dict.get(ID+".5junc.DS.R")]:
            leftHArightPloc= str(feat.location)
        if qualifierDict.get('Sequence') == [RHAprimers_Dict.get(ID+".3junc.DS.F")]:
            rightHAleftPloc= str(feat.location)
        if qualifierDict.get('Sequence') == [RHAprimers_Dict.get(ID+".3gen.R")]:
            rightHArightPloc= str(feat.location)
    #print(leftHAleftPloc,leftHArightPloc,rightHAleftPloc, rightHArightPloc)
    #print(LHAprimers_Dict.get(ID+".5gen.F"))
    #print(LHAprimers_Dict.get(ID+".5junc.DS.R"))
    #print(RHAprimers_Dict.get(ID+".3junc.DS.F"))
    #print(RHAprimers_Dict.get(ID+".3gen.R"))
    if leftHAleftPloc == '':
        return
    
    llstart, llend = leftHAleftPloc.split(":")
    llstart = llstart.strip("[")
    #print(llstart)
    lrstart, lrend = leftHArightPloc.split(":")
    lrend = lrend.strip("] (+)")
    lrend = lrend.strip("] (-)")
    #print(lrend)
    rlstart, rlend = rightHAleftPloc.split(":")
    rlstart = rlstart.strip("[")
    #print(rlstart)
    rrstart, rrend = rightHArightPloc.split(":")
    rrend = rrend.strip("] (+)")
    rrend = rrend.strip("] (-)")

    #print(rrend)   
    LeftAmpliconSeq = genome.seq[int(llstart):int(lrend)]
    RightAmpliconSeq = genome.seq[int(rlstart):int(rrend)]
    Left_Amplicon_GCcontent = GC(LeftAmpliconSeq)
    Right_Amplicon_GCcontent = GC(RightAmpliconSeq)
    print("\nleft amplicon GC content(%): ", round(Left_Amplicon_GCcontent, 4))
    print("\nright amplicon GC content(%): ", round(Right_Amplicon_GCcontent, 4))
    


def isa_group_separator(line):
    return line[0]=='='

class primer_record():
    temp=0

# Function below is unfinished
# Idea was to use one function for both the LHA and RHA primer3 calls
# Would need to set primer_name values based on arm = 'left' or arm = 'right'
# Also set format/output filenames as 'LHA' or 'RHA' accordingly
def primer3call(ID, arm):
    master_found_list=[]
    master_not_found_list=[]
    LHAprimers_0  = {}
    LHAprimers_1  = {}
    LHAprimers_2  = {}
    LHAprimers_3  = {}
    LHAprimers_4  = {}
    LHAprimers_5  = {}
    LHAprimers_6  = {}
    LHAprimers_7  = {}
    LHAprimers_8  = {}
    LHAprimers_9  = {}
    
    

    for try_number in range (1, 5):
        format_file = f"{SCRIPT_DIR}/{ID}_LHA_Primer3_format.txt" #'/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS/'+folder+'/'+ID+'_LHA_Primer3_format.txt' #'Not_found_raw_pass_'+str(counter)+'.txt'    #str(primer3file) #      #This is the primer3 format file for input
        output_file = f'{SCRIPT_DIR}/LHA_Primers_found_round_{try_number}.csv' #'/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS/'+folder+'/LHA_Primers_found_round_'+str(x)+'.csv'           #This file saves all the primers that were not found
        #print("output file: "+ output_file)
        not_found_output = f"{SCRIPT_DIR}/{ID}_LHA_Primer3_output_not_found.txt"
        #print('not found file = {}'.format(not_found_output))
        
        primer_call = f'primer3_core -p3_settings_file primer3_global_settings_{try_number} -output {output_file} {format_file}' # This calls primer3.  You have to load it into the cluster (not possible from VSCode.  Use "module load primer3")
        #print('call to system: {}'.format(primer_call))
        os.system(primer_call)
        
        with open(output_file, 'r') as f:
            for key, group in itertools.groupby(f, isa_group_separator):
                if not key:
                    data={}
                    
                    for item in group:
                        field, value=item.split('=')
                        value=value.strip('\n')
                        data[field]=value
                        #print(item)
                    
                    if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
                        master_found_list.append(data)
                        #print('Master Found List:', master_found_list)
                    else:
                        with open(not_found_output, 'w') as not_found:
                            not_found.write('SEQUENCE_ID='+data['SEQUENCE_ID']+'\n')
                            not_found.write('SEQUENCE_TARGET='+data['SEQUENCE_TARGET']+'\n')
                            not_found.write(str('SEQUENCE_TEMPLATE='+data['SEQUENCE_TEMPLATE']+'\n'))
                            not_found.write('=' + '\n')
                        master_not_found_list.append(data)
                    
                    #print (str(data.get('SEQUENCE_ID'))+"_leftP = "+str(data.get('PRIMER_LEFT_0_SEQUENCE')))
                    #print (str(data.get('SEQUENCE_ID'))+"_rightP = "+str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
                    #print(data.get("SEQUENCE_ID"))
                
                else:
                    pass
        
        #print('Found Primers for {} gRNAs'.format(len(master_found_list)))
        #print('Did not find Primers for {}'.format(len(master_not_found_list)))
        #print(master_found_list)
        #print(master_not_found_list)
        
        found_df=pd.DataFrame(master_found_list)
        #print('Found DF:', found_df)
        found_df['Primer_design_round']=try_number
        found_df.to_csv(output_file)
        
        if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
            #print(str(data['PRIMER_LEFT_NUM_RETURNED']))
            #print(str(data.get('PRIMER_LEFT_0_SEQUENCE')))
            #print(str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
            #for i in range(0, int(data['PRIMER_LEFT_NUM_RETURNED'])):
            
            LHAprimers_0[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_0_SEQUENCE'))
            LHAprimers_0[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_0_SEQUENCE'))
            LHAprimers_1[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_1_SEQUENCE'))
            LHAprimers_1[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_1_SEQUENCE'))
            LHAprimers_2[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_2_SEQUENCE'))
            LHAprimers_2[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_2_SEQUENCE'))
            LHAprimers_3[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_3_SEQUENCE'))
            LHAprimers_3[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_3_SEQUENCE'))
            LHAprimers_4[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_4_SEQUENCE'))
            LHAprimers_4[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_4_SEQUENCE'))

            LHAprimers_5[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_5_SEQUENCE'))
            LHAprimers_5[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_5_SEQUENCE'))
            LHAprimers_6[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_6_SEQUENCE'))
            LHAprimers_6[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_6_SEQUENCE'))
            LHAprimers_7[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_7_SEQUENCE'))
            LHAprimers_7[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_7_SEQUENCE'))
            LHAprimers_8[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_8_SEQUENCE'))
            LHAprimers_8[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_8_SEQUENCE'))
            LHAprimers_9[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_9_SEQUENCE'))
            LHAprimers_9[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_9_SEQUENCE'))

            print("\nLHAprimers:",  LHAprimers_0,'\n',
                                    LHAprimers_1,'\n',
                                    LHAprimers_2,'\n',
                                    LHAprimers_3,'\n',
                                    LHAprimers_4,'\n',
                                    LHAprimers_5,'\n',
                                    LHAprimers_6,'\n',
                                    LHAprimers_7,'\n',
                                    LHAprimers_8,'\n',
                                    LHAprimers_9,'\n')

    #print found_df.head()
    os.remove(not_found_output)
    os.remove(format_file)
    
    return(LHAprimers_0, LHAprimers_1, LHAprimers_2, LHAprimers_3, LHAprimers_4, LHAprimers_5, LHAprimers_6, LHAprimers_7, LHAprimers_8, LHAprimers_9) 

# calls primer3 on LHA primer3 Format file created in homology arm seq
# LHA primer3 Format file exists in homology arm primer finder
def LHAprimer3call(ID, folder):
    #print(primer3file)
    counter = 1
    dict_obj ={}
    master_found_list=[]
    master_not_found_list=[]
    new='\n'
    LHAprimers_0  = {}
    LHAprimers_1  = {}
    LHAprimers_2  = {}
    LHAprimers_3  = {}
    LHAprimers_4  = {}
    LHAprimers_5  = {}
    LHAprimers_6  = {}
    LHAprimers_7  = {}
    LHAprimers_8  = {}
    LHAprimers_9  = {}
    

    for x in range (1, 5):
        if str(ID+".5junc.DS.R") in LHAprimers_4:
            print("stopping LHA primer3call")
            break
        #print(x)
        # master_found_list=[]
        # master_not_found_list=[]
        
        input_file = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_LHA_Primer3_format.txt' #'/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS/'+folder+'/'+ID+'_LHA_Primer3_format.txt' #'Not_found_raw_pass_'+str(counter)+'.txt'    #str(primer3file) #      #This is the primer3 format file for input
        
        #file = open(input_file)
        #print(file.read())
        #file.close()
        
        #LHA_Primers_found = open('/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS/'+folder+'/LHA_Primers_found_round_'+str(x)+'.csv', 'w')
        #LHA_Primers_found.close()
        output_file = 'LHA_Primers_found_round_'+str(x)+'.csv' #'/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS/'+folder+'/LHA_Primers_found_round_'+str(x)+'.csv'           #This file saves all the primers that were not found
        primer_call = 'primer3_core -p3_settings_file primer3_global_settings_'+str(x)+' -output '+output_file + ' ' +input_file          #This calls primer3.  You have to load it into the cluster (not possible from VSCode.  Use "module load primer3")
        #print('call to system: {}'.format(primer_call))
        os.system(primer_call)

        #file = open(output_file)
        #print("output file here")
        #print(file.read())
        #file.close()

        not_found_output = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_LHA_Primer3_output_not_found.txt'
        #print('not found file = {}'.format(not_found_output))
        not_found = open(not_found_output, 'w')
        counter=counter+1
        #print("output file: "+ output_file)
        with open(output_file, 'r') as f:
            for key, group in itertools.groupby(f, isa_group_separator):
                if not key:
                    data={}
                    
                    for item in group:
                        field, value=item.split('=')
                        value=value.strip('\n')
                        data[field]=value
                        #print(item)
                    if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
                        master_found_list.append(data)
                        #print('Master Found List:', master_found_list)
                    
                    else:
                        not_found.write('SEQUENCE_ID='+data['SEQUENCE_ID']+'\n')
                        not_found.write('SEQUENCE_TARGET='+data['SEQUENCE_TARGET']+'\n')
                        not_found.write(str('SEQUENCE_TEMPLATE='+data['SEQUENCE_TEMPLATE']+'\n'))
                        not_found.write('=' + '\n')
                        master_not_found_list.append(data)
                    

                    #print (str(data.get('SEQUENCE_ID'))+"_leftP = "+str(data.get('PRIMER_LEFT_0_SEQUENCE')))
                    #print (str(data.get('SEQUENCE_ID'))+"_rightP = "+str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
                    #print(data.get("SEQUENCE_ID"))
                    
                    
                else:
                    pass

           
        #print('Found Primers for {} gRNAs'.format(len(master_found_list)))
        #print('Did not find Primers for {}'.format(len(master_not_found_list)))
        #print(master_found_list)
        #print(master_not_found_list)
 
        found_df=pd.DataFrame(master_found_list)
                #output_file = str('Primers_'+str(x)+'_pass')
                #print 'output_file= {}'.format(output_file)

        #print('Found DF:', found_df)
        found_df['Primer_design_round']=x
        found_df.to_csv(output_file)
        f.close()
        not_found.close()
        
        
        if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
            #print(str(data['PRIMER_LEFT_NUM_RETURNED']))
            #print(str(data.get('PRIMER_LEFT_0_SEQUENCE')))
            #print(str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
            #for i in range(0, int(data['PRIMER_LEFT_NUM_RETURNED'])):
            
            LHAprimers_0[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_0_SEQUENCE'))
            LHAprimers_0[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_0_SEQUENCE'))
            LHAprimers_1[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_1_SEQUENCE'))
            LHAprimers_1[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_1_SEQUENCE'))
            LHAprimers_2[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_2_SEQUENCE'))
            LHAprimers_2[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_2_SEQUENCE'))
            LHAprimers_3[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_3_SEQUENCE'))
            LHAprimers_3[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_3_SEQUENCE'))
            LHAprimers_4[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_4_SEQUENCE'))
            LHAprimers_4[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_4_SEQUENCE'))

            LHAprimers_5[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_5_SEQUENCE'))
            LHAprimers_5[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_5_SEQUENCE'))
            LHAprimers_6[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_6_SEQUENCE'))
            LHAprimers_6[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_6_SEQUENCE'))
            LHAprimers_7[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_7_SEQUENCE'))
            LHAprimers_7[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_7_SEQUENCE'))
            LHAprimers_8[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_8_SEQUENCE'))
            LHAprimers_8[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_8_SEQUENCE'))
            LHAprimers_9[ID+'.5gen.F'] = str(data.get('PRIMER_LEFT_9_SEQUENCE'))
            LHAprimers_9[ID+'.5junc.DS.R'] = str(data.get('PRIMER_RIGHT_9_SEQUENCE'))

            print("\nLHAprimers:",  LHAprimers_0,'\n',
                                    LHAprimers_1,'\n',
                                    LHAprimers_2,'\n',
                                    LHAprimers_3,'\n',
                                    LHAprimers_4,'\n',
                                    LHAprimers_5,'\n',
                                    LHAprimers_6,'\n',
                                    LHAprimers_7,'\n',
                                    LHAprimers_8,'\n',
                                    LHAprimers_9,'\n')




    #print found_df.head()
    #print("LHA primers here:\n",LHAprimers)
    #for i in range (1, x):
        #source = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/LHA_Primers_found_round_'+str(i)+'.csv'
        #destination = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/'+ID+'_LHA_Primers_found_round_'+str(i)+'_Mark.csv'
        #shutil.copyfile(source, destination)
    os.remove("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_LHA_Primer3_output_not_found.txt')
    os.remove("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_LHA_Primer3_format.txt')
    
    
    return(LHAprimers_0, LHAprimers_1, LHAprimers_2, LHAprimers_3, LHAprimers_4, LHAprimers_5, LHAprimers_6, LHAprimers_7, LHAprimers_8, LHAprimers_9) 




# calls primer3 on RHA primer3 Format file created in homology arm seq
# RHA primer3 Format file exists in homology arm primer finder folder
def RHAprimer3call(ID, folder):
    #print(primer3file)
    counter = 1
    dict_obj ={}
    master_found_list=[]
    master_not_found_list=[]
    new='\n'
    RHAprimers_0 = {}
    RHAprimers_1 = {}
    RHAprimers_2 = {}
    RHAprimers_3 = {}
    RHAprimers_4 = {}
    RHAprimers_5 = {}
    RHAprimers_6 = {}
    RHAprimers_7 = {}
    RHAprimers_8 = {}
    RHAprimers_9 = {}
    

    for x in range (1, 5):
        if str(ID+".3gen.R") in RHAprimers_4:
            print("stopping RHA primer3call")
            break
        #print(x)
        master_found_list=[]
        master_not_found_list=[]
        input_file = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_RHA_Primer3_format.txt' #'Not_found_raw_pass_'+str(counter)+'.txt'    #str(primer3file) #      #This is the primer3 format file for input
        output_file = 'RHA_Primers_found_round_'+str(x)+'.csv'           #This file saves all the primers that were not found
        primer_call = 'primer3_core -p3_settings_file primer3_global_settings_'+str(x)+' -output '+output_file + ' ' +input_file          #This calls primer3.  You have to load it into the cluster (not possible from VSCode.  Use "module load primer3")
        #print('call to system: {}'.format(primer_call))
        os.system(primer_call)
        not_found_output = "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_RHA_Primer3_output_not_found.txt'
        #print('not found file = {}'.format(not_found_output))
        not_found = open(not_found_output, 'w')
        counter=counter+1
        #print("output file: "+ output_file)
        with open(output_file, 'r') as f:
            for key, group in itertools.groupby(f, isa_group_separator):
                if not key:
                    data={}
                    
                    for item in group:
                        field, value=item.split('=')
                        value=value.strip('\n')
                        data[field]=value
                        #print(item)
                    if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
                        master_found_list.append(data)
                    
                    else:
                        not_found.write('SEQUENCE_ID='+data['SEQUENCE_ID']+'\n')
                        not_found.write('SEQUENCE_TARGET='+data['SEQUENCE_TARGET']+'\n')
                        not_found.write(str('SEQUENCE_TEMPLATE='+data['SEQUENCE_TEMPLATE']+'\n'))
                        not_found.write('=' + '\n')
                        master_not_found_list.append(data)
                    
                    #print (str(data.get('SEQUENCE_ID'))+"_leftP = "+str(data.get('PRIMER_LEFT_0_SEQUENCE')))
                    #print (str(data.get('SEQUENCE_ID'))+"_rightP = "+str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
                    #print(data.get("SEQUENCE_ID"))
                else:
                    pass
        #print('Found Primers for {} gRNAs'.format(len(master_found_list)))
        #print('Did not find Primers for {}'.format(len(master_not_found_list)))
        #print(master_found_list)
        #print(master_not_found_list)
 
        found_df=pd.DataFrame(master_found_list)
                #output_file = str('Primers_'+str(x)+'_pass')
                #print 'output_file= {}'.format(output_file)
        found_df['Primer_design_round']=x
        found_df.to_csv(output_file)
        f.close()
        not_found.close()

        if str(data['PRIMER_LEFT_NUM_RETURNED']) >='1':
            
            print(str(data['PRIMER_LEFT_NUM_RETURNED']))
            #print(str(data.get('PRIMER_LEFT_0_SEQUENCE')))
            #print(str(data.get('PRIMER_RIGHT_0_SEQUENCE')))
            #for i in range(0, int(data['PRIMER_LEFT_NUM_RETURNED'])):
            
            RHAprimers_0[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_0_SEQUENCE'))
            RHAprimers_0[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_0_SEQUENCE'))
            RHAprimers_1[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_1_SEQUENCE'))
            RHAprimers_1[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_1_SEQUENCE'))
            RHAprimers_2[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_2_SEQUENCE'))
            RHAprimers_2[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_2_SEQUENCE'))
            RHAprimers_3[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_3_SEQUENCE'))
            RHAprimers_3[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_3_SEQUENCE'))
            RHAprimers_4[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_4_SEQUENCE'))
            RHAprimers_4[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_4_SEQUENCE'))
            
            RHAprimers_5[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_5_SEQUENCE'))
            RHAprimers_5[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_5_SEQUENCE'))
            RHAprimers_6[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_6_SEQUENCE'))
            RHAprimers_6[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_6_SEQUENCE'))
            RHAprimers_7[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_7_SEQUENCE'))
            RHAprimers_7[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_7_SEQUENCE'))
            RHAprimers_8[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_8_SEQUENCE'))
            RHAprimers_8[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_8_SEQUENCE'))
            RHAprimers_9[ID+'.3junc.DS.F'] = str(data.get('PRIMER_LEFT_9_SEQUENCE'))
            RHAprimers_9[ID+'.3gen.R'] = str(data.get('PRIMER_RIGHT_9_SEQUENCE'))

            print("\nRHAprimers:",  RHAprimers_0,'\n',
                                    RHAprimers_1,'\n',
                                    RHAprimers_2,'\n',
                                    RHAprimers_3,'\n',
                                    RHAprimers_4,'\n',
                                    RHAprimers_5,'\n',
                                    RHAprimers_6,'\n',
                                    RHAprimers_7,'\n',
                                    RHAprimers_8,'\n',
                                    RHAprimers_9,'\n'
                                    )
            
    #print found_df.head()
    #print(RHAprimers)
    #for i in range (1, x):
        #source = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/RHA_Primers_found_round_'+str(i)+'.csv'
        #destination = "/research/groups/millergrp/home/common/CORE PROJECTS/"+folder+'/'+ID+'_RHA_Primers_found_round_'+str(i)+'_Mark.csv'
        #shutil.copyfile(source, destination)
    
    os.remove("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_RHA_Primer3_output_not_found.txt')
    os.remove("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_RHA_Primer3_format.txt')
    return(RHAprimers_0, RHAprimers_1, RHAprimers_2, RHAprimers_3, RHAprimers_4, RHAprimers_5, RHAprimers_6, RHAprimers_7, RHAprimers_8, RHAprimers_9) 

        


#uses project ID to find correct folder in CORE PROJECTS, then enumerates possible file options in folder for the user to select from
#returns path to desired file and the name of the folder
def findfile(projectID):
    enumeratedFiles = {}
    enumeratedFolders = {}
    enumerated_filenames = {}
    count = 1
    #directory = "/research/groups/millergrp/home/common/CORE PROJECTS" #"Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS"
    directory = "/research/groups/millergrp/home/common/CORE PROJECTS" #"/research/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/CORE PROJECTS"
    #print(directory)
    #iterate through folders
    for folder in os.listdir(directory):
        #print(folder)
        #find project ID in folder name
        if projectID in folder: #if folder == projectID: #if projectID in folder: 
            #print("Folder Found:",folder,", looking through folder...")
            #print(folder)
            #iterate through files in folder
            for filename in os.listdir(os.path.join(directory, folder)):
                filename_upper = filename.upper()
                #print(filename)
                #check for vector
                if filename_upper.endswith(".GBK") and (filename_upper.__contains__("DONOR") or filename_upper.__contains__("DSDNA") or filename_upper.__contains__("SSDNA") or filename_upper.__contains__("MEGAMER")):
                    #print("found",filename,"in",folder,"folder")
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
            print()
    print()
            
           
    if count == 1:
        print("no files found!!")
        return -1, folder, -1
    else:
        userNum = input("Which of the files above are you using? (type number and return (ex. 2)): ")
    print()
    #print(userNum)
    #print(enumeratedFiles)
    #print(enumeratedFiles.keys())
    #print(enumeratedFiles.get(int(userNum)))
    userfile = enumeratedFiles.get(int(userNum))
    folder_name = enumeratedFolders.get(int(userNum))
    filename1 = enumerated_filenames.get(int(userNum))
    #print("Path to file:     ", userfile)
    return userfile, folder_name, filename1





# change product size range in global settings to correct size 
# these global settings exist in homology arm primer finder

def changeProductSize(HALen):
    #print(HALen)
    print('primer product size ranges: HAlen is: {}, add 10-50, 51-71, then 72-82'.format(HALen))
    for x in range (1, 5):
        with open('primer3_global_settings_'+str(x), 'r+') as file:
            global_line_list = file.readlines()
            global_line_list[0] = 'Primer3 File - http://primer3.sourceforge.net\n'
            global_line_list[3] = 'PRIMER_PRODUCT_SIZE_RANGE= '+str(HALen-85)+'-'+str(HALen)+'\n'
            global_line_list[5] = 'PRIMER_NUM_RETURN=10\n'
            if global_line_list[4] != 'PRIMER_EXPLAIN_FLAG=1\n':
                global_line_list.insert(4,'PRIMER_EXPLAIN_FLAG=1\n')
            file.seek(0)
            file.truncate(0)
            file.writelines(global_line_list)
            file.close()
    file1 = open('primer3_global_settings_1', 'r')
    #print(file1.read())
    file1.close()
    file1 = open('primer3_global_settings_2', 'r')
    #print(file1.read())
    file1.close()
    file1 = open('primer3_global_settings_3', 'r')
    #print(file1.read())
    file1.close()
    file1 = open('primer3_global_settings_4', 'r')
    #print(file1.read())
    file1.close()


# Uses filepath to open genbank file and find homology arms
# Adds 200 bp tails to homology arms 
# Writes the primer3 format to separate files for each homology arm
# Returns RHA and LHA primer3 files and the length of each homology arm
def homologyArmSeq(filepath, ID, folder):
    genome=SeqIO.read(filepath,'genbank') 
    genomeLength = len(genome) - 1
    isLeftHA = False
    isRightHA = False
    
    for feat in genome.features:
        qualifierDict = feat.qualifiers
        if qualifierDict.get('label') == None:
            continue
        if 'LEFT HA' in str(qualifierDict.get('label')).upper() or 'LHA' in str(qualifierDict.get('label')).upper() or str(qualifierDict.get('label')).upper()== ['800BP LHA']:
            print('Left HA location:', feat.location)
            leftHALoc = str(feat.location)
            isLeftHA = True
        if 'Right HA' in str(qualifierDict.get('label')).upper() or 'RHA' in str(qualifierDict.get('label')).upper() or str(qualifierDict.get('label')).upper()== ['800BP RHA']:
            print('Right HA location:', feat.location)
            rightHAloc = str(feat.location)
            isRightHA = True
    if isLeftHA == False or isRightHA == False:
        print("No Homology Arms Found")
        return -1, -1, -1, -1
    LstringList = leftHALoc.split(':')
    LeftHAstart = str(LstringList[0]).replace('[', '')
    LeftHAend = str(LstringList[1]).replace('](+)', '')
    RstringList = rightHAloc.split(':')
    RightHAstart = str(RstringList[0]).replace('[', '')
    RightHAend = str(RstringList[1]).replace('](+)', '')


    Ls = int(LeftHAstart) - 35
    Le = int(LeftHAend) + 50
    Rs = int(RightHAstart) - 50
    Re = int(RightHAend) + 35


    #in case homology arms + the 200 bp tails overlap
    if Le > int(RightHAstart):
        Le = int(RightHAstart) - 1
    if Rs < int(LeftHAend):
        Rs = int(LeftHAend) + 1

    if Ls < 0 or Ls > genomeLength:
        print("not in range")
        return
    if Le < 0 or Le > genomeLength:
        print("not in range") 
        return
    if Rs < 0 or Rs > genomeLength:
        print("not in range")
        return
    if Re < 0 or Re > genomeLength:
        print("not in range")
        return

    #print("Left HA with 200 bp tails: ", Ls, "-", Le)
    #print("Right HA with 200 bp tails: ", Rs, "-", Re)
    print(Ls, Le)
    print(Rs, Re)
    LeftHAseq = genome.seq[Ls:Le]
    RightHAseq = genome.seq[Rs:Re]

    #print(LeftHAseq)
    #print()
    #print(RightHAseq)
    #print()
    
    leftHALen = len(LeftHAseq)
    rightHALen = len(RightHAseq)
    print('LEFT HA is {} bp and is : {}'.format(len(LeftHAseq),LeftHAseq))
    print('RIGHT HA is {}bp and is: {}'.format(len(RightHAseq),RightHAseq))
    #print(rightHALen)

    index = genome.seq.find(RightHAseq, start=0)

    #print("right HA seq index:", index)

    LHAprimer3Format = "SEQUENCE_ID=LeftHAseq\nSEQUENCE_TARGET=34,"+str(leftHALen-73)+"\nSEQUENCE_TEMPLATE="+LeftHAseq+"\n=\n"
    RHAprimer3Format = "SEQUENCE_ID=RightHAseq\nSEQUENCE_TARGET=39,"+str(rightHALen-73)+"\nSEQUENCE_TEMPLATE="+RightHAseq+"\n=\n"
    print('target region for lha: {}'.format(LeftHAseq[34:leftHALen-35]))
    print('target region for rha: {}'.format(RightHAseq[39:rightHALen-35]))
    
    #print('left HA primer3 format:\n', LHAprimer3Format)
    #print('right HA primer3 format:\n', RHAprimer3Format)


    LHAprimer3File = open("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_LHA_Primer3_format.txt', 'w')
    LHAprimer3File.write(str(LHAprimer3Format))
    LHAprimer3File.close()


    RHAprimer3File = open("/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/Donor_Design_Programs/"+ID+'_RHA_Primer3_format.txt', 'w')
    RHAprimer3File.write(str(RHAprimer3Format))
    RHAprimer3File.close()


    return LHAprimer3File, RHAprimer3File, leftHALen, rightHALen



if __name__ == "__main__":
    loop = "Y"
    while True:
        if loop != "Y":
            break
        else:
            jxn_primer_gen()
            loop = input("Would you like to run the Junction Primer Generator again? (y/n): ").upper()
