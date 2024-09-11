import os
import subprocess
import re
import time
import pandas as pd
import argparse


#*Creates a sorted list of guides, guide#s and OTA scores from a CRISPICK output file.
#*End product is passed to library_draft.py to choose how many guides are needed.


parser = argparse.ArgumentParser()

parser.add_argument('-i', help='input design file', default='sgrna-designs.txt', type=str)
parser.add_argument('-l', help='library name', default='test_lib', type=str)
parser.add_argument('-design',help='choose crispick or guidelist to create design',default='crispick',type=str)
parser.add_argument('-s',help='species. (h)uman or (m)ouse',default='h',type=str)
parser.add_argument('-q',help='guide quota/the number of guides per gene',default=5, type=int)
parser.add_argument('-cas',help='cas type choose 9 or 12',default='9',type=str)
parser.add_argument('-ntc',help='type -ntc false if no NTCs are present in crispick design',default=True)
parser.add_argument('-f', help='choose function "md" merge design files, "ml" merge library drafts')


args = parser.parse_args()



cas_type = args.cas
species = args.s
guide_quota = args.q
input_file = args.i
library_name = args.l
need_ntc = args.ntc

draft_file_name = f"{library_name}_draft.xlsx" #*will just be the guides to fill the quota with the best OTA scores
sorted_file_name = f"{library_name}_all_guides.xlsx" #*make sure file is xlsx. Sorted will be master list of all guides and OTA scores 

stars = "*"*60 #for pretty printing


def get_pam(cas_type):
    
    if cas_type == "9":
        pam = 'NNNNNNNNNNNNNNNNNNNNNGG'
    
    elif cas_type == "12":
        pam ='TTTVNNNNNNNNNNNNNNNNNNNNN'
        
    else:
        pam = cas_type

    return pam

def create_offinder_template(input_file, pam, species, crispick_flow):
    
    def combine_columns(rows):
        
        if crispick_flow == True:
            if pam == 'NNNNNNNNNNNNNNNNNNNNNGG':
                if str(rows['Target Gene Symbol']) != '':
                    return str(rows['sgRNA Sequence']) + str(rows['PAM Sequence']) + ' 3 ' + str(rows['Target Gene Symbol'])
                
            elif pam =='TTTVNNNNNNNNNNNNNNNNNNNNN':
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
        input_df = pd.read_csv(input_file, sep='\t')
        input_df['Combined'] = input_df.apply(combine_columns, axis=1)
        input_df = input_df.drop(columns=['Sequence','Name'])

    output_df = pd.concat([pam_df,input_df],ignore_index=True)
    output_df.to_csv("columns_combined_for_offinder.txt",sep='\t',header=False,index=False)    
    
    print("\n\nCasOffinder Template Created\n\n")

def call_offinder():

    input_guides = 'columns_combined_for_offinder.txt'
    
    #load_offinder = f'module load cas-offinder/2.4.1-rhel8'
    offinder_call = f'bsub -P Cas_offinder -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 10000 -o offinder_summary.txt -J cas_offinder_library "hostname & cas-offinder {input_guides} G CasOffinder_Results.txt"'
    
       
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

def create_ota_table(crispick_flow):
    
    positive_controls = ['RPA3','PCNA','DBR1','PLK1','RPL3','KIF11','EEF2','POLR2B','POLR2A','GAPDH','PSMB1']
    
    def _guide_list_capatilization(sorted_df):
        
        #keeps the guide number lower case while allowing for the gene name to be capitalized
        guide_number = sorted_df['name'].str.split(".",n=1,expand=True)
        sorted_df['name'] = guide_number[0]
        sorted_df['guide_num'] = guide_number[1]
        
        #change capatilization of gene name before addding NTC to avoid switchcasing NTC
        if species == 'h':
            sorted_df['name'] = sorted_df['name'].str.upper()
        if species =='m':
            sorted_df['name'] = sorted_df['name'].str.capitalize()
        
        sorted_df['name'] = sorted_df['name'] + "." +sorted_df['guide_num']
        sorted_df = sorted_df.drop('guide_num',axis=1)

        return sorted_df

    def _crispick_capatilization(sorted_df,need_ntc):
        
        def __get_non_target_control():
            
            tmp_df = pd.read_csv(input_file,delimiter="\t",usecols=['Input','sgRNA Sequence'])
            
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

        #change capatilization of gene name before addding NTC to avoid switchcasing NTC
        if species == 'h':
            sorted_df['name'] = sorted_df['name'].str.upper()
        if species =='m':
            sorted_df['name'] = sorted_df['name'].str.capitalize()
        
        #add in NTC's from CRISPICK file, adds NTCs to bototm of list
        
        if need_ntc == True:
            ntc_df = __get_non_target_control()  
            sorted_df = pd.concat([sorted_df,ntc_df])
        else:
            None

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

        arranged_columns = ['Name','gene','full_gRNA','Long_0','Long_1','Long_2','Long_3']
        
        labeled_df = guide_name_df[arranged_columns]
        labeled_df = labeled_df.rename(columns={'name':'Name','gRNA':'Full gRNA'})
        labeled_df = labeled_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})
            
        return labeled_df
    
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

    #crispick designs get the NTCs concated to the end and guide numbers later.  guide_list designs get the guide # split before capitlization
    #seperates the two work flows
    if crispick_flow == True:
        named_df = _crispick_capatilization(sorted_df, need_ntc)
        
    else:
        named_df = _guide_list_capatilization(sorted_df)

    #rearrange and rename columns to fit downstream (library draft) workflow
    arranged_columns = ['name','gRNA','Long_0','Long_1','Long_2','Long_3']
    named_df  = named_df[arranged_columns]
    named_df  = named_df.rename(columns={'name':'Name','gRNA':'full_gRNA'})


    if crispick_flow == True:
        all_guides_df = named_df.pipe(_add_guide_numbers)
    else:
        all_guides_df = named_df
        
    all_guides_df = all_guides_df.astype({"Long_0":'int',"Long_1":'int',"Long_2":'int',"Long_3":'int'})
    
    print(all_guides_df.head())
    
    all_guides_df.to_excel(sorted_file_name,index=False)
    

    print(f"{stars}\nOTA Complete. Sorted OTA all guides list saved as: ' {sorted_file_name} '.\n{stars} ")

def library_drafter(sorted_file_name, guide_quota, draft_file_name):
    
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

    library_df = pd.concat([picked_df,ntc_df],ignore_index=True)
    library_df.index = library_df.index + 1
    
    library_df.to_excel(draft_file_name,index=False)
    
    mageck_df = library_df[['Name','full_gRNA','gene']].copy()
    
    
    mageck_df.to_csv(f'{library_name}_mageck.csv',header=False,index=False)
    
    
    print(stars)
    print(stars)
    print(library_df.head(20))
    print(library_df.tail(20))
    print(stars)
    print(f"Library draft complete.  Please see file {draft_file_name}")
    print(f"Mageck file compiled.  Please see file {library_name}_mageck.csv")
    print(stars)
    
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



def main():
    os.system("module load cas-offinder/2.4.1-rhel8")  
    pam = get_pam(cas_type)
    
    if args.design == 'crispick':
        
        print(f'{stars}')
        print("Generating library design from CRISPick files")
        print(f'{stars}')
        
        create_offinder_template(input_file,pam,species,crispick_flow=True)
        call_offinder()
        create_ota_table(crispick_flow=True)
        library_drafter(sorted_file_name, guide_quota, draft_file_name)
                
        
        return None
        
    elif args.design =='guidelist':
        
        print(f'{stars}')
        print("Generating library design from a list of guides")
        print(f'{stars}')

        create_offinder_template(input_file,pam,species,crispick_flow=False)
        call_offinder()
        create_ota_table(crispick_flow=False)
        #library_drafter(sorted_file_name, guide_quota, draft_file_name)

    clean_files()


main()
