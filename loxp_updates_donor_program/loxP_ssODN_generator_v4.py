from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord
import glob
import os
import re

from add_loxp_to_genbank import add_loxp  # Import function to add loxP to genbank

# Runs w Python3


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
    print(str("{}.anti.ssODN,{}\n".format(name, antisense)))
    print(str("{}.sense.ssODN".format(name)))
    print(str("{}.anti.ssODN\n".format(name)))


def add_primers_to_ssODN_record(ssODN_gene_record, current_gbk, name, project_path):
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
                print("found the reverse complement:{}".format(feature))
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
        ssODN_gene_record.annotations["molecule_type"]="DNA"
        
    #* Fixes issue with no CDSs included in final gbk.  
    
    #generate list of all features to add to loxp record
    all_loxp_features = []
    
    #only take CDS from current gbk
    for feature in current_gbk.features:
        if feature.type == "CDS":
            all_loxp_features.append(feature)
        
    for feature in ssODN_gene_record.features:
        all_loxp_features.append(feature)

    #create loxp SeqIO record from ssODN loxp design
    loxp_seq = Seq(str(ssODN_gene_record.seq))
    loxp_record = SeqRecord(seq=loxp_seq, id="loxp_record.id")
    loxp_record.annotations["molecule_type"]="DNA"
    
    for feature in all_loxp_features:
        loxp_record.features.append(feature)
    
    #calls add_loxp function from add_loxp_to_genbank.py to shift all CDS downstream of Loxp site
    amended_ssODN_gene_record = add_loxp(loxp_record)
    

    SeqIO.write(amended_ssODN_gene_record,str(os.path.join(project_path, f"{name}.loxP.sense.ssODN.gbk")),"genbank",
    )

    print(f"Saved: {name}.loxP.sense.ssODN.gbk")


def find_RE_and_set_loxP_orientation(current_gbk, name, feature, project_path):
    gRNA_orientation_dict = {-1: "bottom", 1: "top"}
    gRNA_cut_site_dict = {
        -1: 6,
        1: 17,
    }  # Position of gRNA cut sites, -1 =nontemplate strand, 1=template strang
    HA_length = 40  # Length of homology arms
    gRNA_cut_site = feature.location.start + gRNA_cut_site_dict[feature.location.strand]
    loxP = "ATAACTTCGTATAGCATACATTATACGAAGTTAT"
    restriction_enzyme_dict = {
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "NotI": "GCGGCCGC",
        "EcoRI": "GAATTC",
    }
    cut_site='None'
    RE_search_region = current_gbk.seq[feature.location.start - 300 : feature.location.start + 330]
    for (enzyme,site) in (restriction_enzyme_dict.items()):  # Check 300bp upstream/downstream for a unique restriction enzyme
        if site not in RE_search_region:
            RE_pick = enzyme
            RE_seq = site
            cut_site=site
            RE_side = str(input("Is this loxP on 5' or 3'? (Enter 5 or 3)->")).strip()
            if RE_side == "5":
                insertion_dict = {
                    "len_1": len(loxP),
                    "label_1": ["loxP","color:31849B"],
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
                    "label_2": ["loxP","color:31849B"],
                }
                loxP_RE = site + loxP
            else:
                print("Not a valid option")
            break
        else:
            pass
    if cut_site=='None':
        print("Warning-no uniqe restriction site found +/-300 bp (enzymes searched: {}".format(str(restriction_enzyme_dict.keys())))

    # Annotate and add features to ssODN GBK file
    full_gene_record_seq = (
        current_gbk.seq[:gRNA_cut_site] + loxP_RE + current_gbk.seq[gRNA_cut_site:]
    )
    full_gene_record = SeqRecord(seq=full_gene_record_seq, id="current_gbk.id")
    full_gene_record.annotations["topology"] = "linear"
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
    display_IDT(name, full_gene_record, gRNA_cut_site, HA_length, len(loxP_RE))
    return full_gene_record


def create_ssODNs_from_input_dir():
    # path = os.path.join(r"Z:\ResearchHome\Groups\millergrp\home\common\Python\CAGE_Programs\LoxP_ssODN_generator\input")
    # save_dir = os.path.normpath(r"Z:\ResearchHome\Groups\millergrp\home\common\Python\CAGE_Programs\LoxP_ssODN_generator\output")

    path = os.path.join(
        r"C:\Users\jconnell\Documents\programming\ngs_design_loxp_ssodn\input"
    )
    save_dir = os.path.normpath(
        r"C:\Users\jconnell\Documents\programming\ngs_design_loxp_ssodn\output"
    )
    os.chdir(path)  # Change path to projects_input folder
    file_list = os.listdir(path)
    print("Creating ssODNs for: {}".format(file_list))
    # Go through each folder in the input dir
    for f in file_list:
        gbks = glob.glob(os.path.join(f, "*mod_NGS.gbk"))
        print("all gbks: {}".format(gbks))
        for mod_NGS in gbks:
            create_loxP_ssODN(mod_NGS, path)


def get_user_input_gRNA(current_gbk):
    # Lists gRNAs in the project mod_NGS.gbk and queries user for choice
    gRNA_dict = {}  # Dict of found gRNAs:features
    for feature in current_gbk.features:
        if feature.type == "misc_feature":
            #pat =  "[A-Z]+\d+[.][A-Z0-9]+([.]g\d+)"  # Pattern to recongize Project_number.gene.gRNA. (i.e. CAGE52.hWASSUP.g1)
            #pat2 = "[A-Z]+\d+[.][A-Z0-9]+[.][A-Z0-9]+[.]+([.]g\d+)"  # Pattern to recongize Project_number.gene.gRNA_type.gRNA. (i.e. CAGE52.hWASSUP.9.g1  or CAGE52.hWASSUP.12a.g1)

            pat =  "[A-Z]+\d+\.[^\.]+\.(g\d+)"  # Pattern to recongize Project_number.gene.gRNA. (i.e. CAGE52.hWASSUP.g1)
            pat2 = "[A-Z]+\d+\.[^\.]+\.([^\.]+)\.(g\d+)"  # Pattern to recongize Project_number.gene.gRNA_type.gRNA. (i.e. CAGE52.hWASSUP.9.g1  or CAGE52.hWASSUP.12a.g1)

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
        gRNA = input("Enter target gRNA # (q to go back) -> ")
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


def create_loxP_from_core_projects_dir():
    # Queries user for project #, checks that only a single project# exists.  Loads the gene_mod_NGS.gbk.
    # Calls get_user_input for gRNA then sends gRNA and file to create loxP
    path = os.path.join(r"/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/")
    #path='\research_jude\\rgs01_jude\groups\millergrp\home\common\CORE PROJECTS'
    os.chdir(path)  # Change path to projects_input folder
    file_list = os.listdir(path)
    project_num = str(input("Enter project # for LoxP ssODN-> ")).upper()
    print("Project#: {}".format(project_num))
    # Go through each folder in the input dir
    for f in file_list:
        #print(f)
        if (project_num == f.split("-")[-1]):  # Check if input project number is equal to a project # in core projects dir
            gbks = glob.glob(os.path.join(f, "*mod_NGS.gbk"))
            path = os.path.join(r"/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/")
            if len(gbks) > 1:
                print("More than one mod_NGS file: {}".format(gbks))
            elif len(gbks) == 0:
                print("No mod_NGS.gbk fils in that project")
            else:
                project_path = path + "/" + f
                mod_NGS = gbks[0]  # Rename gbk to mod_NGS for easier readability
                current_gbk = SeqIO.read(mod_NGS, "genbank")
                gRNA, feature = get_user_input_gRNA(current_gbk)  # query user for target gRNA
                if gRNA != "Q":
                    ssODN_gene_record = find_RE_and_set_loxP_orientation(
                        current_gbk, gRNA, feature, project_path
                    )  # insert loxP and unique restriction enzyme
                    add_primers_to_ssODN_record(
                        ssODN_gene_record, current_gbk, gRNA, project_path
                    )
                    break
                else:
                    break
                # create_loxP_ssODN_from_specified_gRNA(mod_NGS,project_path,path,gRNA)
                # break
            break


def loxP_from_crispr_design(project_dir):
    # This does not work current will probably be abandoned

    path = os.path.join(
        r"/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/"
        + project_dir
    )
    print("loxP generator for: {}".format(project_dir))
    print("path: {}".format(path))
    os.chdir(path)  # Change path to projects_input folder
    file_list = os.listdir(path)
    mod_NGS_gbk = glob.glob("*mod_NGS.gbk")
    print("mod gbks are: {}".format(mod_NGS_gbk))
    for mod_NGS in mod_NGS_gbk:
        create_loxP_ssODN(mod_NGS, path)


def create_loxP_ssODN(mod_NGS, path):
    # path = os.path.join(r"Z:\ResearchHome\ResearchHomeDirs\millergrp\common\Python\CAGE_Programs\LoxP_ssODN_generator\input")
    save_dir = os.path.normpath(
        r"/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS/LoxP_ssODN_generator\output"
    )
    print("Create_loxP_ssODN")
    os.chdir(path)  # Change path to projects_input folder
    current_gbk = SeqIO.read(mod_NGS, "genbank")
    for feature in current_gbk.features:
        if feature.type == "misc_feature":
            print("Feature: {}".format(feature))
            # print(feature)
            pat = "[A-Z]+\d+[.][A-Z]+\d+[.][A-Z0-9]+([.]g\d+)"  # Pattern to recongize Project_number.gene.gRNA. (i.e. CAGE52.hWASSUP.g1)
            try:
                name = str(feature.qualifiers["label"]).strip("['']")
                print("label: {}".format(name))
            except KeyError:
                print("KeyErrorS")
                name = str(feature.qualifiers["note"]).strip("['']")
                print("Note is: {}".format(name))
                # qual_length = len(feature.qualifiers["note"])
                # if qual_length >1:
                #    name = name.split(",")[0].strip("'")
            gRNA_pat = re.match(pat, name)
            if gRNA_pat != None:
                print("Found a gRNA:{}".format(name))
            # qual_length = len(feature.qualifiers["label"])
            ssODN_gene_record = find_RE_and_set_loxP_orientation(
                current_gbk, name, feature
            )
            add_primers_to_ssODN_record(ssODN_gene_record, current_gbk, name)


def main():
    loop = True
    clear = lambda: os.system("cls")
    while loop:
        # print("1: Create_ssODNs_from_input_dir")
        # print("2: Create ssODNs by Project#")
        # clear()
        print("Main menu")
        print("1: Design a loxP ssODN")
        print("Q: quit")
        choice = str(input("Enter choice-> ")).upper()
        if choice == "XX":
            create_ssODNs_from_input_dir()
        elif choice == "1":
            create_loxP_from_core_projects_dir()
        elif choice == "Q":
            loop = False
        else:
            print("Not a valid option")


if __name__ == "__main__":
    main()
