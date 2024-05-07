import pandas as pd
import sys
import os
import glob
import platform
import pathlib as pl
import datetime as dt
import random

#THIS VERSION HAS DATA HARDCODED INTO IT TO MAKE RUNNING FASTER
#Lines 115-127
DEBUGGING = True if len(sys.argv) > 1 and sys.argv[1] in ['-d', 'd', '--d', '--debug', 'debug', '-debug'] else False

class CAGEDirs:
    def __init__(self):
        self.OS_TYPE = platform.system()
        self.CORE_PROJECTS = self.get_CORE_PROJECTS()
        self.NGS = self.get_NGS()
     
    def get_CORE_PROJECTS(self):
        if self.OS_TYPE == 'Linux':
            return '/research_jude/rgs01_jude/groups/millergrp/home/common/CORE\ PROJECTS'
        elif self.OS_TYPE == 'Windows':
            return r'Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS'
        else:
            print('OS_TYPE must be "Linux" or "Windows"...')
            print(f'OS_TYPE is: {self.OS_TYPE}')
            sys.exit(1)
    
    def get_NGS(self):
        if self.OS_TYPE == 'Linux':
            return '/research_jude/rgs01_jude/groups/millergrp/home/common/NGS'
        elif self.OS_TYPE == 'Windows':
            return r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
        else:
            print('OS_TYPE must be "Linux" or "Windows"...')
            print(f'OS_TYPE is: {self.OS_TYPE}')
            sys.exit(1)


class NGSAmplicon:
    def __init__(self, name, proj_loc):
        self.name = name
        self.proj_loc =proj_loc
        self.proj_csv=max(glob.iglob(os.path.join(proj_loc, "*.csv")), key=os.path.getmtime)
        self.start_plate = 0
        self.save_file = os.path.join(proj_loc, "hits.csv")
        self.plate_key = os.path.join(proj_loc, "plate_key.csv")
        self.plate_key_df = pd.DataFrame()
        self.plate_key_dict= {}     #Dictionary that maps NGS plate with the true (1,2,3,4..) plate number for easy comparison
        self.data_df = pd.read_csv(os.path.join(proj_loc, self.proj_csv))
        self.cols_to_add = []
    
    def use_start_plate(self, start_plate):
        self.start_plate = start_plate
        self.data_df['True_key'] = self.data_df['Name'].str.split('-').str[1].str.replace('Plate', '').astype(int) - self.start_plate + 1
        self.data_df['Well'] = self.data_df['Name'].str.split('-').str[2]
        self.data_df['True_key'] = self.data_df.apply(lambda row: f"{row['True_key']}-{row['Well']}", axis = 1)

    def add_selected_columns(self, cols_to_add):
        rename_dict = {}
        for col in cols_to_add:
            rename_dict[col] = f'{self.name}_{col}'
        self.data_df = self.data_df.rename(columns = rename_dict)
        cols_to_add_renamed = [f'{self.name}_{col}' for col in cols_to_add]
        self.cols_to_add = cols_to_add_renamed

def get_projects():
    projects = input(f'Please enter all of the CAGE project names that you would like to include in your search, separated by a space.\n -> ').strip().split()
    return projects

def pick_subdirs_for_project(subdirs, project):
    print(f'{project} subdirs:')
    proj_subdirs = [subdir for subdir in subdirs if project in subdir]
    for idx, subdir in enumerate(proj_subdirs):
        print(f'\t{idx})\t{subdir}')
    
    print(f'Enter the number beside each subdir you would like to include, separated by a space.')
    subdirs_for_extension = [proj_subdirs[int(idx)] for idx in input(' -> ').strip().split()]
    subdirs_str = "\n\t".join(subdirs_for_extension)
    print(f'\nAdding the following subdirs for {project}:\n\t{subdirs_str}')
    print()
    return subdirs_for_extension

def initialize_amplicon(ngs_dir, NGS_date, subdir, DEBUGGING = False, amp_name = '', amp_start = 0):
    subdir_loc = os.path.join(ngs_dir, NGS_date, 'joined', subdir)
    if not DEBUGGING:
        print(f'Regarding {subdir}...')
        amplicon_name = input(f'Short nickname: ').strip()
        
        # Print first and last plates in CSV?
        start_plate = int(input(f'Miller-Plate # to start on: ').strip())
    else:
        amplicon_name = amp_name
        start_plate = amp_start
    
    this_amplicon = NGSAmplicon(amplicon_name, subdir_loc)
    this_amplicon.use_start_plate(start_plate)
    
    idx_TotalIndel = list(this_amplicon.data_df.columns).index('Total_indel')
    potential_columns = this_amplicon.data_df.columns[3:idx_TotalIndel]
    
    if not DEBUGGING:
        print(f'test_list from {amplicon_name}:')
        for idx, col in enumerate(potential_columns):
            print(f'\t{idx})\t{col}')
        print()
        
        print('Would you like to add any of the columns listed above?\nPlease list the number beside each column, separated by a space.')
        cols_to_add_idxs = [int(input_idx) for input_idx in input(' -> ').strip().split()]
        cols_to_add = [potential_columns[idx] for idx in cols_to_add_idxs]
        print(f'\nAdding the following columns: {cols_to_add}\n')
        print()
    else:
        num_cols = min(2, len(potential_columns))
        cols_to_add = random.choices(potential_columns, k = num_cols)
    
    this_amplicon.add_selected_columns(cols_to_add)
    return this_amplicon

def subdirs_to_amplicons(subdirs_to_include, ngs_dir, NGS_date, DEBUGGING = False, names_list = [], start_list = []):
    amplicons_list=[]
    current_idx = 0
    for subdir in subdirs_to_include:
        if not DEBUGGING:
            this_amplicon = initialize_amplicon(ngs_dir, NGS_date, subdir)
        else:
            this_amplicon = initialize_amplicon(ngs_dir, NGS_date, subdir, DEBUGGING, names_list[current_idx], start_list[current_idx])
            current_idx += 1
        amplicons_list.append(this_amplicon)
        print(f'"{this_amplicon.name}" added to plates_list!\n')
    
    return amplicons_list

def create_merged_df(amplicons_list):
    joined_df = pd.DataFrame(columns = ['True_key', 'Total', 'Name'])
    for amplicon in amplicons_list:
        joined_df = joined_df.merge(amplicon.data_df[['True_key', 'Total', 'Name', *amplicon.cols_to_add]], how = 'outer', on = 'True_key', suffixes = [None, f'_{amplicon.name}'])
    joined_df = joined_df.drop(columns = ['Total', 'Name'])
    col_ordering = ['True_key']
    col_ordering.extend([col for col in joined_df.columns if 'Total' in col])
    col_ordering.extend([col for col in joined_df.columns if 'Total' not in col and 'Name' not in col and col != 'True_key'])
    col_ordering.extend([col for col in joined_df.columns if 'Name' in col])
    return joined_df, col_ordering

def save_merged_df(amplicons_list, joined_df, col_ordering):
    def get_col_widths(dataframe):
        # First we find the maximum length of the index column
        # idx_max = max([len(str(s)) for s in dataframe.index.values] + [len(str(dataframe.index.name))])
        # Then, we concatenate this to the max of the lengths of column name and its values for each column, left to right
        return [
            max([len(str(s)) for s in dataframe[col].values] + [len(col) + 2])
            for col in dataframe.columns
        ]
    
    if not DEBUGGING:
        print('Choose an output folder:')
        for idx, amplicon in enumerate(amplicons_list):
            print(f'{idx})\t{amplicon.name}\n\t{pl.Path(amplicon.proj_loc).stem}')
        print()
        output_dir = amplicons_list[int(input("Number -> ").strip())].proj_loc
        output_name = input("Output Name -> ").strip() + '.xlsx'
        print(f'Outputting to: \n\t{os.path.join(output_dir, output_name)}')
        print()
    else:
        output_dir = amplicons_list[0].proj_loc
        out_time = dt.datetime.now().strftime('%y%m%d-%H%M%S')
        output_name = f'debugging-{out_time}.xlsx'
    
    with pd.ExcelWriter(os.path.join(output_dir, output_name)) as writer:
        joined_df[col_ordering].to_excel(
            writer, columns=col_ordering, sheet_name="hits", index=False,
        )
        worksheet = writer.sheets["hits"]
        for i, width in enumerate(get_col_widths(joined_df[col_ordering])):
            worksheet.set_column(i, i, width)
    print(f'File saved as:\n\t{os.path.join(output_dir, output_name)}\n')
    return


def main():
    cdirs = CAGEDirs()
    ngs_dir = cdirs.NGS
    
    if not DEBUGGING:
        NGS_date = input("Enter the NGS date -> ").strip()
        projects = get_projects()
        
        subdirs = next(os.walk(os.path.join(ngs_dir, NGS_date, 'joined')))[1] # Get all subdir/project names from this NGS dir
        subdirs_to_include = []
        for project in projects:
            picked_subdirs = pick_subdirs_for_project(subdirs, project)
            subdirs_to_include.extend(picked_subdirs)
        print(subdirs_to_include)
        
        amplicons_list = subdirs_to_amplicons(subdirs_to_include, ngs_dir, NGS_date)
    else:
        print(f'Entering debug mode (preset input) - Argument provided: {sys.argv[1]}')
        NGS_date = "012423"
        names_list=['SPAST', 'ABE','ssODN']
        subdirs_to_include=['CAGE1045_hSPAST_F_R_FLAG_ssODN','CAGE719_hPPP1R12C_F_R_short_ABE','CAGE1127_hTANGO1_F_R_test_ssODN']
        start_list=[9,49,17]
        projects = ['CAGE1045', 'CAGE719', 'CAGE1127']
        
        amplicons_list = subdirs_to_amplicons(subdirs_to_include, ngs_dir, NGS_date, DEBUGGING, names_list, start_list)

    joined_df, col_ordering = create_merged_df(amplicons_list)
    
    save_merged_df(amplicons_list, joined_df, col_ordering)
    
    print("Comparison complete")


if __name__ == "__main__":
    main()
