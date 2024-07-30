import os
import subprocess
import re
import time
import pandas as pd


#*Creates a sorted list of guides, guide#s and OTA scores from a CRISPICK output file.
#*End product is passed to library_draft.py to choose how many guides are needed.


input_file = "lib140-sgrna-designs.txt"

draft_file_name = "140_draft.xlsx" #*will just be the guides to fill the quota with the best OTA scores
sorted_file_name = "140_sorted.xlsx" #*make sure file is xlsx. Sorted will be master list of all guides and OTA scores 

cas_type = '9'
species = 'h' #h = human, m = mouse.  Can expand in future if need be.

guide_quota = 5

#* Don't touch these
custom_type = None

#TODO Don't forget to change this if using a short list of guides that do not include ntc's
need_ntc = True


positive_controls = ['RPA3','PCNA','DBR1','PLK1','RPL3','KIF11','EEF2','POLR2B','POLR2A','GAPDH','PSMB1']
stars = "*"*60

def get_pam(cas_type):
    
    if cas_type == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif cas_type == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = cas_type

    return pam

def combine_columns(rows):
    
    if pam == 'NNNNNNNNNNNNNNNNNNNNNGG':
        if str(rows['Target Gene Symbol']) != '':
            return str(rows['sgRNA Sequence']) + str(rows['PAM Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])
        
    elif pam =='TTTVNNNNNNNNNNNNNNNNNNNNN':
        if str(rows['Target Gene Symbol']) != '':
            return str(rows['PAM Sequence']) + str(rows['sgRNA Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])

def format_crispick(input_file, pam, species):
    
    if species.lower() == 'm':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/mm10'
    
    if species.lower() == 'h':
        genome = '/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/bowtie_indexes/fasta/hg38'
        
    pam_df = pd.DataFrame([[genome],[pam]],columns=['Combined'])

    input_df = pd.read_csv(input_file, sep='\t', usecols=['Target Gene Symbol','sgRNA Sequence','PAM Sequence'])
    
    input_df['Combined'] = input_df.apply(combine_columns, axis=1)

    input_df = input_df.drop(columns=['Target Gene Symbol','sgRNA Sequence'])
    input_df.drop(columns=['PAM Sequence'], inplace=True)

    output_df = pd.concat([pam_df,input_df],ignore_index=True)

    output_df.to_csv("columns_combined_for_offinder.txt",sep='\t',header=False,index=False)    
        
    #make seperate df for NTCs

def offinder():

    input_guides = 'columns_combined_for_offinder.txt'
    
    offinder_call = f'bsub -P Cas_offinder -notify done -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 10000 -o run_summary.txt -J cas_offinder "hostname & cas-offinder {input_guides} G CasOffinder_Results.txt"'
                    
    Long_sout = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')    
    jid_Long9 = re.split("<|>", Long_sout)[1]
    if jid_Long9.isnumeric():
        print(stars)
        print("Starting CasOffinder\n")
        print("System will notify when the job is complete")
    else:
        print("Failed in submitting job for casOffinder")
        
    os.system(f"bwait -w 'ended(cas_offinder)'")
    print("CasOffinder Complete.  Formatting OTA table")

def create_ota_table():
    casoff_file = 'CasOffinder_Results.txt'

    #summary frame serves as the template for the completed mismatch dataframe
    tmp_df =pd.read_csv(casoff_file,delimiter="\t", names=["gRNA","chrom", "pos", "result", "strand", "mismatch", "name"], header=None, index_col=False)

    tmp_df.drop(['chrom', "strand", 'result'],axis=1,inplace=True)
    
    #get gene names now to merge later
    gene_df = tmp_df[['gRNA','name']]
    
    #groups df by the position, mismatch and gRNA seq.  Counts the number of mismatches per position and creates
    #multi-index df that lists the positions as columns for counting across the series        
    count_df = tmp_df[['gRNA','mismatch','pos']].groupby(by=['gRNA', 'mismatch']).count().unstack()

    #dropping one of the indexes for a flatter df
    count_df.columns = count_df.columns.droplevel(1)
    count_df.fillna(0, inplace=True)
    
    count_df.columns = ['Long_0','Long_1','Long_2','Long_3']
    
    #counts across the series of positional mismatches and renames them Long_0 etc
    count_df['Long_3'] = count_df.iloc[:,0:4].sum(axis=1)
    count_df['Long_2'] = count_df.iloc[:,0:3].sum(axis=1)
    count_df['Long_1'] = count_df.iloc[:,0:2].sum(axis=1)
    
    #reset index to get rid of any multi-index shannigans from above step and merge with gene names
    #gets rid of duplicatese.  count_df will not have any NTC's.  
    count_df.reset_index(inplace=True)
    count_df = count_df.merge(gene_df, on="gRNA", how="right",right_index=False).drop_duplicates()
    #ignore index to reset it
    sorted_df = count_df[['gRNA','name','Long_0','Long_1','Long_2','Long_3']].sort_values(by=['name','Long_0','Long_1','Long_2','Long_3'],ignore_index=True)

    #change capatilization of gene name before addding NTC to avoid switchcasing NTC
    if species == 'h':
        sorted_df['name'] = sorted_df['name'].str.upper()
    if species =='m':
        sorted_df['name'] = sorted_df['name'].str.capitalize()
    
    #add in NTC's from CRISPICK file, adds NTCs to bototm of list
    
    if need_ntc == True:
        ntc_df = get_non_target_control()
        #sorted_df.sort_values(by=['Long_0'],ascending=True,inplace=True)    
        sorted_df = pd.concat([sorted_df,ntc_df])
    else:
        None
    
    #rearrange and rename columns to fit downstream (library draft) workflow
    arranged_columns = ['name','gRNA','Long_0','Long_1','Long_2','Long_3']
    sorted_df = sorted_df[arranged_columns]
    sorted_df = sorted_df.rename(columns={'name':'Name','gRNA':'full_gRNA'})

    sorted_df = (sorted_df.pipe(add_guide_numbers))
    
    print(sorted_df.head())
    
    sorted_df.to_excel(sorted_file_name,index=False)
    

    print(f"{stars}\nOTA Complete. Sorted OTA List saved as: ' {sorted_file_name} '.\n{stars} ")
    
def get_non_target_control():
    ntc_file = input_file
    tmp_df = pd.read_csv(ntc_file,delimiter="\t",usecols=['Input','sgRNA Sequence'])
    
    #select all negative control sequences and rename them NTC
    ntc_df = tmp_df.loc[tmp_df['Input'] == '(NEG_CONTROL)']
    ntc_df = ntc_df.rename(columns={'Input': 'name','sgRNA Sequence': 'gRNA'})
    
    #append pam onto end of NTC sequence
    ntc_df['gRNA'] = ntc_df['gRNA'].astype(str) + 'NGG'
    
    #ntc_df = ntc_df.loc[ntc_df['gRNA'].str.contains('NTC')]
    
    ntc_df.reset_index(inplace=True)
    ntc_df.drop(['index'],axis=1,inplace=True)
    
    #change name of (NEG_CONTROL) to NTC
    ntc_df = ntc_df.replace('(NEG_CONTROL)','NTC')
    #add in Long_0 etc
    ntc_df[['Long_0','Long_1','Long_2','Long_3']] = pd.DataFrame([[0,0,0,0]], index=ntc_df.index)

    return ntc_df

def add_guide_numbers(sorted_df):
    input_df = sorted_df
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

    arranged_columns = ['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3']
    
    labeled_df = guide_name_df[arranged_columns]
    labeled_df = labeled_df.rename(columns={'name':'Name','gRNA':'Full gRNA'})
    labeled_df = labeled_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})
    
    return labeled_df

    #print(lib_df.head(10))


    #lib_df.to_excel(sorted_file_name,header=True,index=False)

def library_drafter(sorted_file_name, guide_quota, ntc_df,output_file):
      
    #read in xlsx as df
    sorted_df = pd.read_excel(sorted_file_name)
    
    
    picked_list=[]
    ntc_list=[]
    cur_gene=""
    prev_gene=""
    guide_count=0
    
    for row in sorted_df.itertuples():
        cur_gene=row[2]
        
        if prev_gene == cur_gene and cur_gene != "NTC" and guide_count < guide_quota:#adding until quota filled
            picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            prev_gene = cur_gene
            guide_count +=1
            
        elif prev_gene == cur_gene and cur_gene != "NTC" and guide_count >= guide_quota: #quota hit
            None
            
        elif cur_gene == "NTC":#parse out NTC df
            ntc_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            
        else: #first time new gene encountered
            guide_count=1
            picked_list.append((row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
            prev_gene=cur_gene
            
    
    picked_df = pd.DataFrame(picked_list,columns=['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3'])
    ntc_df = pd.DataFrame(ntc_list,columns=['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3'])
    picked_df.index = picked_df.index+1
    ntc_df.index = ntc_df.index+1
    
    

    picked_df.sort_values(by=['Name','Long_0','Long_1','Long_2','Long_3'],inplace=True)
    #Creates a series from Name column (gene.g#).  Creates lambda that splits the seires @ '.g', takes the #side and converts that series as int.
    #s stays a series until the very end when its cast as in with astype.  Have to treat the entire thing as a series
    if need_ntc == True:
        try:
            ntc_df.sort_values(by=['Name'],inplace=True, key=lambda s: (s.str.rsplit('.g',expand=True)[1]).astype(int))
        except KeyError as key_error:
            print(key_error)
            print("Make sure the input list has NTC's picked.")
    else:
        None
    #print("Library Draft: ")
    #print(picked_df.tail(20))

    library_df = pd.concat([picked_df,ntc_df],ignore_index=True)
    library_df.index = library_df.index + 1
    
    
    print(stars)
    print(stars)
    print(library_df.head(20))
    print(library_df.tail(20))
    print(stars)
    print(f"Library Draft Complete.  Please see file {output_file}")
    print(stars)
    
    
    library_df.to_excel(output_file,index=False)
    
    return library_df.shape[0]

def clean_files():
    #waits to make sure all files are done writing before going to next step
    time.sleep(3)
    
    trash_files = [
        "run_summary.txt",
        "columns_combined_for_offinder.txt",
        "pysub.log",
        "CasOffinder_Results.txt"
    ]
    for file in trash_files:
        try:
            os.remove(file)
        except:
            None



    
if custom_type != None:
    cas_type = custom_type

pam = get_pam(cas_type)
ntc_df = get_non_target_control()

format_crispick(input_file,pam,species)

offinder()

create_ota_table()

#get_non_target_control()

library_drafter(sorted_file_name, guide_quota,ntc_df,draft_file_name)

clean_files()
