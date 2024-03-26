import os
import pandas as pd
import shutil

inputCSV = "input.csv"
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
holding_dir = os.path.join(os.getcwd(),'holding')


all_indel_list=[]

def get_projects():

    clear_holding_dir()
    
    project_list=[]   
    
    #get input from csv
    input_df = pd.read_csv(inputCSV)

    input_projects = input_df.values.tolist()
    #print(project_list)

    for proj in input_projects:
    
        cage_proj, ngs_date, plate_num = proj
        
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined')
        os.chdir(joined_dir)
        
        target_py = cage_proj+'.py'
        
        shutil.copy(target_py,str(holding_dir))
        
        target = [cage_proj,plate_num]
        
        project_list.append(target)
        
    return project_list

def rewrite_target_py(project_list):
    test_list_open = False
    os.chdir(holding_dir)
    
    the_Seq =""
    the_Seq_start=""
    the_Seq_end=""
    fastq_files=""
    test_list=""
    
    for project in project_list:
        
        script, plate_num = project
        
        fastq_files = plate_range(plate_num)
        
        #TODO
        '''
        
        
        with open(script+'.py',"r+",encoding='utf-8') as f:
            code = f.readlines()
        f.close()
        
        for line in code:
            line = line.strip()

            #Find seq, seq_start, seq_end, test
            #find the plate num, redo directory, and change to clone check main

            if line.startswith('the_Seq '):
                the_Seq = line

            elif line.startswith('the_Seq_start'):
                the_Seq_start = line

            elif line.startswith('the_Seq_end'):
                the_Seq_end = line
                
            elif line.startswith('fastq_files'):
                fastq_files = line
                
            elif line.startswith('test_list'):
                test_list_open = True
                test_list = line
                
            elif line.startswith("(") and test_list_open == True:
                test_list = test_list + line
                
            elif line.startswith("]"):
                test_list = test_list + line
                test_list_open = False

            template=f'''

#TODO
    '''

import os,sys,inspect
import argparse

def main():
    ID = str(os.path.basename(__file__)).split('.py')[0]
    {the_Seq}
    {the_Seq_start}
    {the_Seq_end}
    fastq_files = '{fastq_files}'
    {test_list}
    
    search_fastq(ID,the_Seq,the_Seq_start,the_Seq_end,fastq_files,test_list)
    
    
if __name__ =='__main__':		# use this if you want to include modules from a subfolder
        cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0], r'Z:\\ResearchHome\\Groups\\millergrp\\home\\common\\Python\\CAGE_Programs\\NGS_one_main_program')))
        if cmd_subfolder not in sys.path:
            sys.path.insert(0, cmd_subfolder)
            from clone_check_main_all_indels import search_fastq
        main()
    '''
    #TODO
    '''
                
        new_py = str(script).replace(".py","") + "_all_indels.py"

        with open(new_py,'w') as f:
            for line in template:
                f.write(line)
        f.close()
    
        all_indel_list.append(new_py)

    for indel_script in all_indel_list:
        run_updated_scripts(indel_script)
    
    #TODO
    '''



def run_updated_scripts(indel_script):
    os.chdir()
    
    os.system(f'python {indel_script}')
    
    return


def clear_holding_dir():
        
    trash_files = os.scandir(holding_dir)
    for file in trash_files:
        os.remove(file)    
    
def plate_range(plate_num):
    
    plate_range_expression = None
    
    #read in plate numbers
    
    plates = list(plate_num.split(","))
    
    print(f"The plates {plates}")
    
    i=0
    if len(plates) > 1:
        base_len = len(plates[0])
        if len(set(map(len,plates)))!=base_len:
            print("Plate matched")
            print(f"here are the plates: {plates}")
        else:
            print(f"plates different digits {plates}")
    else:
        print(f"Single plate: {plates}")
    #seperate them
    #figure out if how to run on plates 19-20?
    
    
    
    
    

    
    return plate_range_expression
#parse project details from csv

project_list = get_projects()


rewrite_target_py(project_list)









#TODO how to check for multiple plates and create a plate range to run 