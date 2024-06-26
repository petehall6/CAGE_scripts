import os
import pandas as pd
import shutil

inputCSV = "input.csv"
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'
holding_dir = os.path.join(os.getcwd(),'holding')
all_indel_list=[]

def find_and_run_projects():

    clear_holding_dir()
    
    target_list=[]   
    
    #get input from csv
    input_df = pd.read_csv(inputCSV)

    #puts columns into list so that I can containerize the data eaiser
    input_projects = input_df.values.tolist()

    for proj in input_projects:
        
        cage_proj, ngs_date, plate_num = proj
        
        joined_dir = os.path.join(ngs_dir,ngs_date,'joined')
        os.chdir(joined_dir)
        
        target_py = cage_proj+'.py'
        
        shutil.copy(target_py,str(holding_dir))
        
        target = [cage_proj,plate_num]
        
        target_list.append(target)
        
        rewrite_target_py(target_list,joined_dir)
        #empties the target_list after every rewrite
        #again this is for containerization
        target_list.pop(0)

def rewrite_target_py(project_list,joined_dir):
    
    test_list_open = False
    os.chdir(holding_dir)
    
    the_Seq =""
    the_Seq_start=""
    the_Seq_end=""
    fastq_files=""
    test_list=""

    
    for project in project_list:
        
        script, plate_num = project
        fastq_files,multi_range = plate_range(plate_num)
        
        #reads in the old file and parses out the seq, seq_start, seq_end and corrects the NGS_one_main directory if need be.  
        with open(script+'.py',"r+",encoding='utf-8') as f:
            code = f.readlines()
        f.close()
        
        for line in code:
            line = line.strip()

            if line.startswith('the_Seq '):
                the_Seq = line

            elif line.startswith('the_Seq_start'):
                the_Seq_start = line

            elif line.startswith('the_Seq_end'):
                the_Seq_end = line
            
            #if multiple guides are present in the original code, this keeps the brackets open while it fills in the list
            elif line.startswith('test_list'):
                test_list_open = True
                test_list = line
                
            elif line.startswith("(") and test_list_open == True:
                test_list = test_list + line
                
            elif line.startswith("]"):
                test_list = test_list + line
                test_list_open = False
            #template for just one plate/range
            single_template=f'''
            
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
    
            #template for the multi_range clone check main
            multi_template=f'''
#TODO make multi_range template


import os,sys,inspect
import argparse

def main():
    ID = str(os.path.basename(__file__)).split('.py')[0]
    {the_Seq}
    {the_Seq_start}
    {the_Seq_end}
    #fastq_files = *.fastq
    {test_list}
    
    
    glob_range = {fastq_files}
    
    search_fastq(ID,the_Seq,the_Seq_start,the_Seq_end,glob_range,test_list)
    
    
if __name__ =='__main__':		# use this if you want to include modules from a subfolder
        cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0], r'Z:\\ResearchHome\\Groups\\millergrp\\home\\common\\Python\\CAGE_Programs\\NGS_one_main_program')))
        if cmd_subfolder not in sys.path:
            sys.path.insert(0, cmd_subfolder)
            from clone_check_main_multi_all_indels import search_fastq
        main()
    '''
        

        new_py = str(script).replace(".py","") + "_all_indels.py"

        if multi_range == False:
            with open(new_py,'w') as f:
                for line in single_template:
                    f.write(line)
            f.close()
        
        else:
            with open(new_py,'w') as f:
                for line in multi_template:
                    f.write(line)
            f.close()

        run_updated_scripts(new_py, joined_dir)
        
def run_updated_scripts(indel_script,joined_dir):
    
    os.chdir(holding_dir)
    shutil.copy(indel_script,joined_dir)
    
    os.chdir(joined_dir)
    
    print(f"\n\nRunning: {indel_script}\n\n")
    
    os.system(f'python {indel_script}')
    
def clear_holding_dir():
    trash_files = os.scandir(holding_dir)
    for file in trash_files:
        os.remove(file)    

def get_ones_place(plate_list):
    #helper function for plate_range
    
    for plate in plate_list:

        ones = plate[-1]
            
    return ones

def plate_range(plate_num):
    
    #read in and sanitize plate numbers
    plates = list(str(plate_num).strip().replace(" ","").split(","))
    
    #adds 0 before single digit numbere
    for index in range(0, len(plates)):
        if len(plates[index]) < 2:
            plates[index] = f'0{plates[index]}'
    
    tens_match=False
    multi_range=False
    if len(plates) > 1:
        base_len = len(plates[0])
        #creates a map of the lengths of plates numbers for each entry of csv.
        #Set removes duplicates so there's only 1 value to compare
        #seperates 00-99 from 100-1XX
        
        if len(set(map(len,plates))) != base_len:
            #print(f"Plate same digits: {plates}")
            #same tens in two digits
            if base_len == 2:
                for index in range(1, len(plates)):
                    if plates[0][0] == plates[index][0]:
                        tens_match = True
                        #continue to ones place parsing below                        
                    else:
                        #run multi range
                       # print(f"Different tens?: {plates[0][0]} : {plates[index][0]}")
                        tens_match = False
                        #first mismatch breaks to run multi range
                        break
                        
                if tens_match == True:
                    #gets the ones place for all the plates in plate list and sorts them
                    #maps are really cool
                    ones = sorted((map(get_ones_place, plates)))
                    updated_fastq_files = f"Miller-Plate{plates[0][0]}[{ones[0]}-{ones[-1]}]*.fastq"
                    print(updated_fastq_files)
                        
                else:
                    #print("tens didnt match in two digit")
                    multi_range = True
                    glob_range=[]
                    
                    for plate in plates:
                        glob_range.append(f"Miller-Plate{plate}*.fastq")

                    updated_fastq_files = glob_range
                    print(updated_fastq_files)

                    

            #same tens in three digits
            if base_len == 3:
                print("triple digits")
                for index in range(1, len(plates)):
                    if plates[0][1] == plates[index][1]:
                        tens_match = True
                        #continue to ones place parsing below                        
                    else:
                        #run multi range
                        print(f"Different tens: {plates[0][1]} : {plates[index][1]}")
                        tens_match = False
                        #first mismatch breaks to run multi range
                        break
                        
                if tens_match == True:
                    #gets the ones place for all the plates in plate list and sorts them
                    #maps are really cool
                    ones = sorted((map(get_ones_place, plates)))
                    #parse first two digits of 3 digit number '[:2]
                    updated_fastq_files = f"Miller-Plate{plates[0][:2]}[{ones[0]}-{ones[-1]}]*.fastq"
                    print(updated_fastq_files)
                        
                else:
                    print("tens didnt match in triple digit")
                    multi_range = True
                    glob_range=[]
                    
                    for plate in plates:
                        glob_range.append(f"Miller-Plate{plate}*.fastq")

                    updated_fastq_files = glob_range
                    print(updated_fastq_files)

        else:
            print("Run multi_range")
            multi_range = True
            glob_range=[]

            for plate in plates:
                glob_range.append(f"Miller-Plate{plate}*.fastq")

            updated_fastq_files = glob_range
            print(updated_fastq_files)
        
        
    else:
        print(set(map(len,plates)))
        updated_fastq_files = f"Miller-Plate{plates[0]}*.fastq"
        print(f"Single plate: {plates}")
        print(updated_fastq_files)
        multi_range=False
        
    
    return updated_fastq_files, multi_range

def main():
    find_and_run_projects()
    clear_holding_dir()
    
main()