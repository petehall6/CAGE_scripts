import os
import subprocess
import re
import time
import pandas as pd
import argparse
from pathlib import Path
import shutil
from tabulate import tabulate
from add_cloning_extensions_offinder import add_extension


#*Creates a sorted list of guides, guide#s and OTA scores from a CRISPICK output file.
#*End product is passed to library_draft.py to choose how many guides are needed.
stars = "*"*60 #for pretty printing


parser = argparse.ArgumentParser()

parser.add_argument('-i','--input', help='input design file.  Can be crispick file or a list of guides.  Include file extensions.', default='sgrna-designs.txt', type=str)
parser.add_argument('-l','--lib', help='library name', default='libxxx', type=str)
parser.add_argument('-d','--design', help='choose crispick or list to create design',type=str)
parser.add_argument('-s','--species', help='species. (h)uman or (m)ouse', default='h',type=str)
parser.add_argument('-q','--quota', help='guide quota/the number of guides per gene',default=5, type=int)
parser.add_argument('-c','--cas', help='cas type choose 9 or 12',default='9', type=str)
parser.add_argument('-n','--ntc',help='Enter percentage or exact number of NTCs. Add "%%" after if using percentage. Default value of 10%%.  enter 0 if you dont want any'.format(),default='10%')
parser.add_argument('-f','--func',help='choose function "merge_crispick" merge sgrna-design files, \n"update_mageck" update mageck file from draft, \n"twist" add twist primers, "merge_drafts" merge library drafts',default='None')
parser.add_argument('--cwd', help='print working directory', type=str)

args = parser.parse_args()

print(stars)
for arg, value in vars(args).items():
    print(f"{arg}: {value}")
print(stars+"\n\n")

design_type = args.design
cas_in = str(args.cas)
species_in = args.species
guide_quota = args.quota
input_file = args.input
library_name = args.lib
ntc_percentage = args.ntc
aux_func = args.func
bash_dir = args.cwd

SCRIPT_DIR = Path(__file__).resolve().parent
MERGE_DIR = Path(bash_dir).joinpath('merge_folder')


draft_file_name = f"{bash_dir}/{library_name}_draft.xlsx" #*will just be the guides to fill the quota with the best OTA scores
sorted_file_name = f"{bash_dir}/{library_name}_all_guides.xlsx" #*make sure file is xlsx. Sorted will be master list of all guides and OTA scores 

def get_species(species):
    
    if species.lower() == 'mouse':
        return 'm'
    if species.lower() == 'human':
        return 'h'
    
    elif species.lower() == 'm' or species.lower() == 'h':
        return species.lower()
    else:
        print("Species not recognized.  Please enter 'mouse', 'm','human', or 'h'")
        exit()
    
def get_cas_type(cas_type):
    
    if cas_type.lower() == 'cas9':
        return '9'
    if cas_type.lower() == 'cas12':
        return '12'
    
    elif cas_type.lower() == '9' or cas_type.lower() == '12':
        return cas_type.lower()
    else:
        print("Cas not recognized.  Please enter 'cas9', '9', 'cas12', '12'")
        exit()

def get_pam(cas_type):
    
    if cas_type == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif cas_type == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = cas_type

    return pam

def create_offinder_template(input_file, pam, species, crispick_flow):
    
    def combine_columns(rows):
        
        if crispick_flow == True:
            if pam == 'NNNNNNNNNNNNNNNNNNNNNGG':
                if str(rows['Target Gene Symbol']) != '':
                    return str(rows['sgRNA Sequence']) + str(rows['PAM Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])
                
            elif pam =='TTTVNNNNNNNNNNNNNNNNNNNNNNN':
                if str(rows['Target Gene Symbol']) != '':
                    return str(rows['PAM Sequence']) + str(rows['sgRNA Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])
        else:
            return rows['Sequence']+ ' 3 ' +rows['Name']
            
    if species.lower() == 'm':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/mm10'
    
    if species.lower() == 'h':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/hg38'
        
    pam_df = pd.DataFrame([[genome],[pam]],columns=['Combined'])

    if crispick_flow == True:
        input_df = pd.read_csv(input_file, sep='\t', usecols=['Target Gene Symbol','sgRNA Sequence','PAM Sequence'])
        input_df['Combined'] = input_df.apply(combine_columns, axis=1)
        input_df = input_df.drop(columns=['Target Gene Symbol','sgRNA Sequence'])
        input_df.drop(columns=['PAM Sequence'], inplace=True)
    else:
        input_df = pd.read_excel(input_file, engine='openpyxl')
        print(input_df)
        input_df['Combined'] = input_df.apply(combine_columns, axis=1)
        input_df = input_df.drop(columns=['Sequence','Name'])

    output_df = pd.concat([pam_df,input_df],ignore_index=True)
    
    output_df.to_csv("accessory_files/columns_combined_for_offinder.txt",sep='\t',header=False,index=False)    
    
    print("\n\nCasOffinder Template Created\n\n")

def call_offinder():

    input_guides = 'accessory_files/columns_combined_for_offinder.txt'
    
    #load_offinder = f'module load cas-offinder/2.4.1-rhel8'
    offinder_call = f'bsub -P Cas_offinder -q gpu_priority -gpu "num=1/host" -R a100 -M 10000 -o offinder_summary.txt -e offinder_error_log.txt -J cas_offinder_library "hostname & cas-offinder {input_guides} G accessory_files/CasOffinder_Results.txt"'

    send_job = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')    
    job_id = re.split("<|>", send_job)[1]
    if job_id.isnumeric():
        print(stars)
        print("Starting CasOffinder\n")
        print("System will notify when the job is complete")
    else:
        print("Failed in submitting job for casOffinder")
        
    os.system(f"bwait -w 'ended(cas_offinder_library)'")
    print("CasOffinder Complete.  Formatting OTA table")

def create_ota_table(crispick_flow,species,cas_type):
    
    positive_controls = ['RPA3','PCNA','DBR1','PLK1','RPL3','KIF11','EEF2','POLR2B','POLR2A','GAPDH','PSMB1']
    
    def _guide_list_capatilization(sorted_df):
        
        
        def __guide_list_ntcs(sorted_df):
            
            try:
                all_ntc_df = pd.read_csv(input_file)
            except:
                all_ntc_df = pd.read_excel(input_file,engine='openpyxl')
            
            all_ntc_df = all_ntc_df[all_ntc_df['Name'].str.contains('NTC')]
            
            offinder_ntc_df = sorted_df[['name','gRNA','Long_0','Long_1','Long_2','Long_3']]
            
            offinder_ntc_df = offinder_ntc_df.loc[offinder_ntc_df['name'].str.contains('NTC')]
            
            offinder_ntc_df.reset_index(inplace=True)
            offinder_ntc_df.drop(['index'],axis=1,inplace=True)
            offinder_ntc_df.rename(columns={'gRNA':'Sequence','name':'Name'},inplace=True)
            #offinder_ntc_df[['Long_0','Long_1','Long_2','Long_3']] = pd.DataFrame([[0,0,0,0]], index=offinder_ntc_df.index)
            
            
            #*lots of name swapping but I dont have time to fix it.  I'll just keep it as is
            ntc_df = pd.concat([offinder_ntc_df,all_ntc_df]).drop_duplicates(subset=['Sequence']).reset_index(drop=True)
            ntc_df.fillna(0,inplace=True)
            ntc_df.rename(columns={'Sequence':'gRNA','Name':'name'},inplace=True)
            
            return ntc_df
        
        
        ntc_df = __guide_list_ntcs(sorted_df)

        sorted_df = pd.concat([sorted_df,ntc_df]).drop_duplicates(subset=['gRNA']).reset_index(drop=True)
        
        gene_list = sorted_df['name'].values.tolist()

        #append positive control tag to positive controls based on name
        for gene in gene_list:
            if gene in positive_controls:
                gene_list[gene_list.index(gene)] = str(gene + "_pos")

        guide_name_list = []
        cur_gene =''
        prev_gene = ''
        g_num = 1

        for gene in gene_list:
            cur_gene = gene
            if cur_gene == prev_gene:
                guide_name_list.append(str(gene) + '.g' + str(g_num))
                prev_gene = gene
                g_num +=1
            else:
                g_num = 1
                guide_name_list.append(str(gene) + '.g' + str(g_num))
                prev_gene = gene
                g_num +=1

        #replaces Name with Gene.g#
        guide_name_df = sorted_df.assign(Name=guide_name_list)

        arranged_columns = ['gRNA','name','Name','Long_0','Long_1','Long_2','Long_3']
        
        labeled_df = guide_name_df[arranged_columns]
        labeled_df = labeled_df.rename(columns={'name':'Gene',"Name":'name'})
        labeled_df = labeled_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})

        #keeps the guide number lower case while allowing for the gene name to be capitalized
        guide_number = labeled_df['name'].str.split(".",n=1,expand=True)
        labeled_df['name'] = guide_number[0]
        labeled_df['guide_num'] = guide_number[1]
        
        #change capatilization of gene name before addding NTC to avoid switchcasing NTC
        if species == 'h':
            labeled_df['name'] = labeled_df['name'].str.upper()
        if species =='m':
            labeled_df['name'] = labeled_df['name'].str.capitalize()
        
        labeled_df['name'] = labeled_df['name'] + "." +labeled_df['guide_num']
        labeled_df = labeled_df.drop('guide_num',axis=1)
        
        return labeled_df

    def _crispick_capatilization(sorted_df):

        #change capatilization of gene name before addding NTC to avoid switchcasing NTC
        if species == 'h':
            sorted_df['name'] = sorted_df['name'].str.upper()
        if species =='m':
            sorted_df['name'] = sorted_df['name'].str.capitalize()

        return sorted_df
    
    def _add_guide_numbers(named_df):
        input_df = named_df
        pos_count = 0
        input_df['gene'] = input_df['Name']
        gene_list = input_df['Name'].values.tolist()

        #append positive control tag to positive controls based on name
        for gene in gene_list:
            if gene in positive_controls:
                gene_list[gene_list.index(gene)] = str(gene + "_pos")

        #print(gene_list)
        guide_name_list = []
        cur_gene =''
        prev_gene = ''
        g_num = 1

        for gene in gene_list:
            cur_gene = gene
            if cur_gene == prev_gene:
                guide_name_list.append(str(gene) + '.g' + str(g_num))
                prev_gene = gene
                g_num +=1
            else:
                g_num = 1
                guide_name_list.append(str(gene) + '.g' + str(g_num))
                prev_gene = gene
                g_num +=1

        #replaces Name with Gene.g#
        guide_name_df = input_df.assign(Name=guide_name_list)

        arranged_columns = ['Name','gene','gRNA','Long_0','Long_1','Long_2','Long_3']
        
        labeled_df = guide_name_df[arranged_columns]
        labeled_df = labeled_df.rename(columns={'name':'Name'})
        labeled_df = labeled_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})
            
        return labeled_df
    
    casoff_file = 'accessory_files/CasOffinder_Results.txt'

    #summary frame serves as the template for the completed mismatch dataframe
    tmp_df =pd.read_csv(casoff_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None, index_col=False)

    tmp_df.drop(['chrom', "strand", 'result'],axis=1,inplace=True)
    #get gene names now to merge later
    gene_df = tmp_df[['gRNA','name']]
    
    #groups df by the position, mismatch and gRNA seq.  Counts the number of mismatches per position and creates
    #multi-index df that lists the positions as columns for counting across the series        
    count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()


    count_df.columns = count_df.columns.droplevel(0)
    count_df.reset_index(0)
    count_df.fillna(0, inplace=True)
    
    #checks for missing columns and fills in the blank
    columns_needed = set([0,1,2,3])
    column_set = set(count_df.columns)
    missing_cols = list(sorted(columns_needed - column_set))
    for column in missing_cols:
        count_df.insert(loc=column, column=column, value=[0 for i in range(count_df.shape[0])])
    
    count_df.columns = ['Long_0','Long_1','Long_2','Long_3']
    #add Long_3 first so as not to change the value of Long_0,Long_1 and Long_2 before adding up to Long_3
    count_df['Long_3'] = count_df.iloc[:,0:4].sum(axis=1)
    count_df['Long_2'] = count_df.iloc[:,0:3].sum(axis=1)
    count_df['Long_1'] = count_df.iloc[:,0:2].sum(axis=1)
    count_df['Long_0'] = count_df.iloc[:,0:1].sum(axis=1)

    #reset index to get rid of any multi-index shannigans from above step and merge with gene names
    #gets rid of duplicatese.  count_df will not have any NTC's.  
    count_df.reset_index(inplace=True)
    count_df = count_df.merge(gene_df, on="gRNA", how="right",right_index=False).drop_duplicates()
    #ignore index to reset it
    sorted_df = count_df[['gRNA','name','Long_0','Long_1','Long_2','Long_3']].sort_values(by=['name','Long_0','Long_1','Long_2','Long_3'],ignore_index=True)


    if cas_type == '9':
        sorted_df['gRNA'] = sorted_df['gRNA'].str.slice(0,20)
    elif cas_type == '12':
        sorted_df['gRNA'] = sorted_df['gRNA'].str.slice(4,27)

    #crispick designs get the NTCs concated to the end and guide numbers later.  guide_list designs get the guide # split before capitlization
    #seperates the two work flows
    if crispick_flow == True:
        named_df = _crispick_capatilization(sorted_df)
        
    else:
        #Offinder will not return data for NTC's that have all 0's.  Will only return hits.
        #Need to parse all NTC's from input list and see what's missing from sorted_df
        #Append 0's scores to NTC's that are missing from sorted_df

        named_df = _guide_list_capatilization(sorted_df)

    #rearrange and rename columns to fit downstream (library draft) workflow
    arranged_columns = ['name','gRNA','Long_0','Long_1','Long_2','Long_3']
    
    named_df  = named_df[arranged_columns]
    
    named_df  = named_df.rename(columns={'name':'Name','gRNA':'gRNA'})
    
    if crispick_flow == True:
        all_guides_df = named_df.pipe(_add_guide_numbers)
    else:
        named_df['gene'] = named_df['Name'].str.replace(r'\.g\d+','',regex=True)
        all_guides_df = named_df[['Name','gene','gRNA','Long_0','Long_1','Long_2','Long_3']]
        
    all_guides_df = all_guides_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})
    
    print(all_guides_df)
    
    all_guides_df.to_excel(sorted_file_name,index=False)
    

    print(f"{stars}\nOTA Complete. Sorted OTA all guides list saved as: ' {sorted_file_name.split('/')[-1]} '.\n{stars} ")

def library_drafter(sorted_file_name, species, cas_type, guide_quota, ntc_percentage, draft_file_name):

    def _get_ntc_df(lib_size, cas_type, species, ntc_percentage):                

        if cas_type == '9':
        
            if species == 'h':
                tmp_df = pd.read_csv(r'accessory_files/ntc_9_human_master_list.txt',delimiter="\t")
                ntc_number_max = 2774 
                
            elif species == 'm':
                tmp_df = pd.read_csv(r'accessory_files/ntc_9_mouse_master_list.txt',delimiter="\t")
                ntc_number_max = 2322
                
        elif cas_type == '12':
            
            if species == 'h':
                tmp_df = pd.read_csv(r'accessory_files/ntc_12_human_master_list.txt',delimiter="\t")
                ntc_number_max = 2000
            elif species == 'm':
                tmp_df = pd.read_csv(r'accessory_files/ntc_12_mouse_master_list.txt',delimiter="\t")
                ntc_number_max = 1999
        
        #parse correct number of ntcs
        #number of validated human ntc is 2775
        #arg can be either whole number or percentage is passed as string
        if  "%" in ntc_percentage:
            ntc_amount = int(ntc_percentage.strip('%')) / 100
            ntc_number = int(round((lib_size * ntc_amount),0)) 
        else:
            ntc_number = int(ntc_percentage)
        
        if ntc_number > ntc_number_max:
            ntc_number = ntc_number_max
        else:
            ntc_number = int(ntc_number)
        
        print(f"Number of NTCs: {ntc_number}")
        
        
        ntc_df = tmp_df.head(ntc_number)
        
        print(f"NTC dataframe: {ntc_df}")

        return ntc_df
    
    #read in xlsx as df
    sorted_df = pd.read_excel(sorted_file_name)
    
    picked_list=[]
    cur_gene=""
    prev_gene=""
    guide_count=0
    
    if guide_quota != 0:
        for row in sorted_df.itertuples():
            cur_gene=row[2]
            
            if prev_gene == cur_gene and guide_count < guide_quota:#adding until quota filled
                picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
                prev_gene = cur_gene
                guide_count +=1
                
            elif prev_gene == cur_gene and guide_count >= guide_quota: #quota hit
                None
                
            else: #first time new gene encountered
                guide_count=1
                picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
                prev_gene=cur_gene
    #Will catch guide list designs since they dont have a quota
    else:
        picked_list = sorted_df.values.tolist()    

    picked_df = pd.DataFrame(picked_list,columns=['Name','gene','gRNA','Long_0','Long_1','Long_2','Long_3'])
    picked_df.index = picked_df.index+1
    picked_df.sort_values(by=['Name','Long_0','Long_1','Long_2','Long_3'],inplace=True)
        
    lib_size = len(picked_df)
    
    ntc_df = _get_ntc_df(lib_size, cas_type, species, ntc_percentage)

    library_df = pd.concat([picked_df,ntc_df],ignore_index=True)
    library_df.index = library_df.index + 1
    
    library_df.fillna(0,inplace=True)
    
    #drops PAM
    if cas_type == '9': 
        library_df['gRNA'] = library_df['gRNA'].str.slice(0,20)
        
    elif cas_type == '12':
        library_df['gRNA'] = library_df['gRNA'].str.slice(4,27)
    
    library_df.to_excel(draft_file_name,index=False)
    
    counts = library_df.value_counts('gene')
    counts.to_excel(f'{library_name}_guide_counts.xlsx',header=False)
    
    mageck_df = library_df[['Name','gRNA','gene']].copy()
    mageck_df.to_csv(f'{library_name}_mageck.csv',header=False,index=False)
    
    if library_df.shape[0] < 20:
        guides_to_print = 5
    else:
        guides_to_print = 20
    
    
    print(stars)
    print(stars)
    print(library_df.head(guides_to_print))
    print(library_df.tail(guides_to_print))
    print(stars)
    print(f"Library draft complete.  Please see file {str(draft_file_name).split('/')[-1]}\n")
    print(f"Mageck file compiled.  Please see file {library_name}_mageck.csv\n")
    print(stars)
    
    return library_df.shape[0]

def merge_crispick_files():
    
    os.chdir(MERGE_DIR)

    crispick_df_list=[]

    #get all files in sgrna-designs.txt
    for file in Path().iterdir():
            #print(file.name)
            if file.name != 'merged_crispick_designs.txt':
                tmp_df = pd.read_csv(file,delimiter="\t")
                print(f'{file.name}:number of rows {tmp_df.shape[0]}')
                crispick_df_list.append(tmp_df)
    
    print(f"Found {len(crispick_df_list)} files to merge")
    print(f"Creating merged_crispick_designs.txt")

    combined_df = pd.concat(crispick_df_list, ignore_index=True)
    
    combined_df.to_csv(f'{library_name}_merged_crispick_designs.txt',sep='\t',index=False)
    
    print(combined_df)
    
    #trying out pathlib to move files around
    source = Path(f'{library_name}_merged_crispick_designs.txt')
    dest = Path(bash_dir).joinpath(f'{library_name}_merged_crispick_designs.txt')
    
    shutil.copy(source,dest)
    
    os.chdir(bash_dir)
    return None

def merge_drafts():

    library_name = input("Library Name: ")
    drafts = input("Please leave off the file extenesions. \nEnter draft files names separated by a comma: ")
    
    input_files = drafts.replace(" ","").split(",")
        
    combined_df = pd.DataFrame()
    
    for file in input_files:
        tmp_df = pd.read_excel(file+".xlsx",engine='openpyxl')
        combined_df = pd.concat([combined_df,tmp_df])

    final_df = combined_df.drop_duplicates()
    
    print(final_df)
    
    final_df.to_excel(f'{library_name}_merged_drafts.xlsx',index=False)
    
    print(stars)
    print(f"Drafts merged.  Merged draft saved as {library_name}_merged_drafts.xlsx")
    print(stars)
    
    update_mageck()
    
    return None

def update_mageck():
    os.chdir(bash_dir)
    draft_file = input("Enter draft file name: ").strip()
    
    mageck_name = input("Enter file name for the new updated mageck file name: ").strip()
    
    draft_df = pd.read_excel(draft_file,engine='openpyxl')
    
    mageck_df = draft_df[['Name','gRNA','gene']].copy()
    
    mageck_df.to_excel(mageck_name+".xlsx",index=False,header=False)
    print(stars)
    print("Update complete")
    print(stars)
    return None

def add_twist_primers():
    
    vector_confirmation = ''
    
    lib_name = input("Enter library name: ").strip().lower()
    
    #check to see if draft file exists
    while draft_exists == False:
        input_filename = input("\n\nEnter the filename of the library draft that will be used to create the oligo pool, include the file extension: ").strip().lower()
        draft_exists = os.path.exists(input_filename)
    
    
    #reads column B from lib_vectors_and_primers.xlsx
    vector_primer_excel = '/research_jude/rgs01_jude/groups/millergrp/home/common/Screens/lib_vectors_and_primers.xlsx'
    
    backbone_df = pd.read_excel(vector_primer_excel, sheet_name='extension_for_TWIST', engine='openpyxl', usecols="A,B,C,D")
    
    backbone_df = backbone_df[backbone_df['cloning_strategy'].str.contains('gibson')].reset_index(drop=True)
    backbone_df.index = backbone_df.index + 1
    
    
    while vector_confirmation != 'y' and vector_confirmation != 'yes':
        print("\n\n")
        print(tabulate(backbone_df, headers='keys',tablefmt='grid',showindex=True))
        print("\n")
            
        line_choice = int(input("Please enter vector line number: "))-1
        
        vector_choice = backbone_df.iloc[line_choice]['Vector']
        
        print(f"Vector Choice: {vector_choice}")
        
        vector_confirmation = input(f"Confirm vector choice: {vector_choice} (y/n): ").strip().lower()
    
    
    desired_primer = backbone_df.iloc[line_choice]['primer_set_name']
    
    input(f"Desired primers are: {desired_primer}")
    
    print(f"Adding {desired_primer} primers to {lib_name}")
    
    add_extension(lib_name, desired_primer, input_filename)
    
    return None

def clean_files():
#waits to make sure all files are done writing before going to next step
    time.sleep(1)
    trash_files = [
        "offinder_summary.txt",
        "columns_combined_for_offinder.txt",
        "offinder_error_log.txt",
        "pysub.log",
        "CasOffinder_Results.txt"
    ]
    os.chdir(bash_dir)
    for file in trash_files:
        try:
            os.remove(file)
        except:
            None
    os.chdir(SCRIPT_DIR)
    for file in trash_files:
        try:
            os.remove(file)
        except:
            None

def main(design_type, cas_in, species_in, guide_quota, input_file, ntc_percentage, aux_func): 
    
    cas_type = get_cas_type(cas_in)
    species = get_species(species_in)
    pam = get_pam(cas_type)
    
    
    if aux_func == 'merge_crispick':
        print(f'{stars}')
        print("Merging crispick files")
        merge_crispick_files()
        print(f'{stars}')
        
        return None
    
    if aux_func == "update_mageck":
        print(f'{stars}')
        print("Updating mageck file")
        print(f'{stars}')
        
        update_mageck()
        
        return None
    
    if aux_func == "twist":
        print(f'{stars}')
        print("Adding TWIST Primers")
        print(f'{stars}')
        
        add_twist_primers()
        
        return None
    
    if aux_func == "merge_drafts":
        print(f'{stars}')
        print("Merging drafts")
        print(f'{stars}')
        
        merge_drafts()
        
        return None
    
    if design_type == 'crispick':
        print(f'{stars}')
        print("Generating library design from CRISPick files")
        print(f'{stars}')
        
        crispick_flow=True
        
        create_offinder_template(input_file,pam,species,crispick_flow)
        call_offinder()
        create_ota_table(crispick_flow, species,cas_type)
        library_drafter(sorted_file_name, species, cas_type, guide_quota, ntc_percentage, draft_file_name)

    if design_type =='list':
        print(f'{stars}')
        print("Generating library design from a list of guides")
        print(f'{stars}')

        crispick_flow=False
        guide_quota = 0

        create_offinder_template(input_file,pam,species,crispick_flow)
        call_offinder()
        create_ota_table(crispick_flow, species,cas_type)
        library_drafter(sorted_file_name, species, cas_type, guide_quota, ntc_percentage, draft_file_name)
    
    clean_files()


main(design_type, cas_in, species_in, guide_quota, input_file, ntc_percentage, aux_func)