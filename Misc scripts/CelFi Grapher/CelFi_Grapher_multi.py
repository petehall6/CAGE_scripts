#7/6/2023 PMH

#This script plots Fitness Assay NGS Data

import os
import pandas as pd
import shutil
import matplotlib as mpl
import matplotlib.pyplot as plt
import importlib as ilib
import subprocess
import glob




#*STEP ONE
#*If you need to run all-indels
#* MOVE PYTHON SCRIPTS TO cell_fitness/python_scripts folder
#*make changes to plate numbers so it runs faster



NGS_DATE = input("Enter NGS run date: ")
#NGS_DATE = '082223'
NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"
working_dir = NGS_DIR+"/"+NGS_DATE+"/joined/"

if os.path.exists(working_dir+'/_cell_fitness_folder') == False:
    os.chdir(working_dir)
    print("Making Folders")
    os.mkdir("_cell_fitness_folder")
    os.mkdir("_cell_fitness_folder/python_scripts")
    os.mkdir("_cell_fitness_folder/csv_files")
    os.mkdir("_cell_fitness_folder/csv_files/ratios")
    
def _run_all_indels():
    #to run or not to run?
    run_all_indels = input("Do you want to run all_indels program?  y or n: ").upper()
    #run_all_indels = 'Y'
    if run_all_indels == 'Y':
        #find project folder
        
        target_dirs=[]
        #TODO make a version that loops through csv to run all stuff.
        
        #***DROP THE PYTHON SCRIPTS (SHORTS,F_R,WHATEVER) IN THE _cell_fitness_folder
        #*Script will create all_indels_
        
        #ngs_date = input("Enter NGS run date: ")
        print(working_dir)
        try: 
            os.chdir(working_dir)
            print("Date found")
        except:
            print("Date not found")
                
        #gets input csv with list of folders
        fit_folder = working_dir+'_cell_fitness_folder/python_scripts'
        
        os.chdir(fit_folder)
        #add all_indels_no_round to clone_check_main
        target_scripts = []
        for f in os.scandir(fit_folder):
            if f.name.endswith(".py"):
                target_scripts.append(f.name)
                
                
        for script in target_scripts:
            os.chdir(fit_folder)
            new_script = script.replace(".py",'')+"_all_indels_no_round.py"
            
            all_indels_script = str(shutil.copy(script,new_script))


            code_replace=[]
            with open(all_indels_script, "rt") as script:
                code = script.readlines()
                for line in code:
                    if line.startswith("    fastq_files = "):
                        #TODO figure out a way to narrow down to specific plates
                        line = line.replace(line, "    fastq_files = 'Miller-Plate*.fastq'\n" )
                    if "from clone_check_main" in line:
                        line = line.replace("clone_check_main", "clone_check_main_all_indels_no_round")
                    code_replace.append(line)

            with open(all_indels_script, 'w+') as script:
                script.writelines(code_replace)        
                
            #copy to joined folder and run
            new_script = shutil.copyfile(all_indels_script, working_dir+'/'+all_indels_script)
            target_dirs.append(new_script)
            os.chdir(working_dir)
            print("Running all_indels.py now.")
            subprocess.run(['python',new_script])
            
            csv_folder = working_dir+"_cell_fitness_folder/python_scripts"
            csv_list = []
            for file in target_dirs:
                os.chdir(working_dir)
                script_dir = file.replace(".py","/")
                os.chdir(script_dir)
                for file in os.scandir(script_dir):
                    if file.name.endswith('all_indels.csv'):
                        csv_list.append(file.name)
                        
                for file in csv_list:
                    print("Copying to csv_folder")
                    shutil.copy(file,csv_folder+file)
            input("All indel programs have finished running.  Check .CSV files and move to the CSV folder.  Press Enter to continue")
    else:
        None


def _input(working_dir, NGS_DATE):
    #get input
    #run_date = "062023"
    run_date = NGS_DATE


    #check if run_date is correct/folder exsists
    while os.path.isdir(working_dir) == True:
        print("Success!  Run date found.")
        break
        
    else:
        print("Run date not found.  Please reenter run date.")
        run_date = input("Enter run date: ").strip()
        working_dir = os.path.join(NGS_DIR, run_date, "joined").replace("/","\\")
        
    
    cagenum_dirs =[]
    #check project number/cage number.  If there anything in the cagenum/projectnum dirs list loop breaks and contiues to emailer
    while len(cagenum_dirs) == 0:
        #cage_num = '1111'
        cage_num = input("Please enter project number: ").upper().strip()
        print("Searching for project folder.  This may take a moment.")
        cagenum_dirs = [dir.name for dir in os.scandir(working_dir) if dir.is_dir() if cage_num in dir.name]
        print(f"Number of projects found: {len(cagenum_dirs)}")
        #empty list evaluate to false
        if not cagenum_dirs:
            print(f"Did not find project number {cage_num}.")
                    
    return working_dir,cagenum_dirs,cage_num

def _rename_csv(working_dir, cagenum_dirs):
    
    csv_dir = os.path.join(working_dir,"_cell_fitness_folder/csv_files/")
    
    
    #create a df from cagenum_dirs for easier reading
    temp_df = pd.DataFrame(cagenum_dirs, columns=['Directory'])
    #set first row to display as 1
    dirs_df = temp_df
    temp_df.drop(columns=temp_df.columns)
    dirs_df.index +=1
    print()
    print(dirs_df)
    print()
    #index still initilazes to 0. Subtract 1 to account for print/actual format
    dir_choice = int(input("Enter directory choice number: "))-1
    #dir_choice = 0
    print()
    cage_dir = dirs_df.iloc[dir_choice,0]

    print(f"CAGE Directory: {cage_dir}")

    #point back to NGS folder to scrape fastq file names from CSV
    target_dir = os.path.join(working_dir,cage_dir)
    os.chdir(target_dir)
    #for saving bargraphs via interactive window
    graph_dir = str(target_dir) + '\\graphs'
    
    if os.path.exists(graph_dir) == False:
        os.mkdir(graph_dir)
    mpl.rcParams["savefig.directory"] = graph_dir

    #finds all .csv's in target dir by looking at file extension. Doesn't need to be a list but it's easier to work with.
    csv_list = [ f.name for f in os.scandir(target_dir) if f.name.endswith(".csv") ]
    temp_df = pd.DataFrame(csv_list, columns=['File'])
    csv_list_df = temp_df
    del temp_df
    csv_list_df.index +=1
    print()
    print(csv_list_df)
    print()
    
    csv_choice = int(input("Please select CSV: "))-1
    #csv_choice = 0
    target_csv = os.path.join(os.path.join(target_dir,csv_list[csv_choice]))
    
    input(f"target_csv: {target_csv}")
    
    #reads CSV into a df to parse plate names
    raw_df = pd.read_csv(target_csv)
    
    csv_df = raw_df[['Sample','0bp','In-frame','Out-of-frame']]
    csv_df.index +=1
    print()
    print(csv_df)
    
    rename = input("Do you wish to rename samples? (y)es or (n)o:  ").upper()
    
    
    if rename == 'Y' or rename == 'YES':
    
        #rename rows.  print vertical and add backk to csv_df
        sample_list =[]
    
        index=1
        #get sample names
        while len(sample_list) < csv_df.shape[0]:       
            sample_name = input(f"Please enter sample name for row {index}: ")
            sample_list.append(sample_name)
            index = index + 1

        
        #create temp_df of sample names for easier merging
        temp_df = pd.DataFrame(sample_list,columns=['Sample'])
        raw_df['Sample'] = temp_df['Sample']
        del temp_df
        #set index to name stop addtional column from being added
        raw_df.set_index(['Name'], inplace=True)
            
        #write sample names to csv
        raw_df.to_csv(os.path.join(target_dir,csv_list[csv_choice]))
        input(f"saving to {csv_dir}")
        raw_df.to_csv(os.path.join(csv_dir,csv_list[csv_choice]))
        del raw_df
        
        print("\nNaming complete\n\n")
        print(target_csv)
        return(target_csv,graph_dir)
    else:
        raw_df.to_csv(os.path.join(target_dir,csv_list[csv_choice]))
        input(f"saving to {csv_dir}")
        raw_df.to_csv(os.path.join(csv_dir,csv_list[csv_choice]))
        print(target_csv)
        return(target_csv,graph_dir)

def _make_plot():
    
    def label_bar(bars):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2.0, 1.0*height,
                            f'{height}',
                            ha='center', va='bottom')
    
    csv_folder = working_dir+"/_cell_fitness_folder/csv_files"
    
    print("hello I should be making a plot")
    
    target_files = []
    for f in os.scandir(csv_folder):
        if f.name.endswith('.csv'):
            target_files.append(f)
    
    
    print(target_files)
    for target_csv in target_files:
    
    
        guide_name = []
        ratio = []
        
        named_csv = pd.read_csv(target_csv)
        oof_df = named_csv[['Sample','0bp','In-frame','Out-of-frame']]
        
        oof_df.index +=1
        print(oof_df)
        print()
        print()
        
        
        
        gene_name = str(target_csv).split("_")[1]

        gene_correct = input(f"Confirm gene name is correct {gene_name}, y or n: ").upper()
        
        print(f"GENE CORRECT: {gene_correct}")
        
        if gene_correct == 'N' or gene_correct == 'NO':
            gene_name = input("Please enter gene name: ")

        
        #gene_name="TESTGENE"
        chart_title = gene_name +" Cell Fitness (CelFi) Assay"
        
        #num_of_guides = 1
        num_of_guides = int(input("How many guides are being compared (1 or 2): "))
        
        if num_of_guides == 2:
            #prompt for time poitns

            
            #TODO. FOR TESTING
            #guide_name = ['Guide One','Guide Two']
            #firstguide_first_time = 2
            #firstguide_last_time = 3
            #secondguide_first_time = 4
            #secondguide_last_time = 5
            
            
            firstguide_first_time = int(input("\nEnter row number of first time point for first guide: "))
            firstguide_last_time = int(input("\nEnter row number of last time point for first guide: "))
            
            secondguide_first_time = int(input("\nEnter row number of first time point for second guide: "))
            secondguide_last_time = int(input("\nEnter row number of last time point for second guide: "))
            
            
            guide_name.append(str(oof_df.at[firstguide_first_time,'Sample']).split(' ')[0])
            guide_name.append(str(oof_df.at[secondguide_first_time,'Sample']).split(' ')[0])
            
            
            print(f"1st guide: {guide_name[0]}")
            print(f"2nd guide: {guide_name[1]}")
            
            guide_name_good = input("Please confirm guide names, y or n: ").upper()
            
            if guide_name_good == 'N' or guide_name_good == 'NO':
                guide_name.clear()
                guide_name.append(input("1st guide name: "))
                guide_name.append(input("2nd guide name: "))
            
            ratio.append(round(oof_df.at[firstguide_last_time, 'Out-of-frame']/oof_df.at[firstguide_first_time, 'Out-of-frame'],2))
            
            ratio.append(round(oof_df.at[secondguide_last_time, 'Out-of-frame']/oof_df.at[secondguide_first_time, 'Out-of-frame'],2))
            
                
            print()
            print()
            print(f"{guide_name[0]} fitness ratio: {ratio[0]}")
            print()
            print()
            print(f"{guide_name[1]} fitness ratio: {ratio[1]}")
            
            #sets overall graph as subplot for easier editing
            fig, ax = plt.subplots()
            #makes the individual guide bar plot a separate entity for easier editing
            
            
            
            #Okay, so this is a super ugly hacky way of doing this....
            #x_pos sets the column number for each bar plot.  In the case of two guides there are two empty columns flanking the guides on either side [ blank, g1, g2, blank]
            #this allows for the bars to be more centered in the chart rather than super far apart
            x_pos=[0,1,2,3]
            bars = ax.bar(x_pos,
                        # sets values for empty columns as 0 and explicitly calls guide score indices
                        [0,ratio[0],ratio[1],0], 
                        color='#16A085', 
                        edgecolor='#17202A',
                        width=0.5,
                        )


            label_bar(bars)
            #uses a short range of 1,2 to get two middle x_pos columns and uses guide_name list to fill in names
            plt.xticks(range(1,3), guide_name)
            
            #lazy way of compiling ratios for csv export
            g1_ratio = [gene_name, guide_name[0], ratio[0]]
            g2_ratio = [gene_name, guide_name[1], ratio[1]]
            
            ratio_list = [g1_ratio,g2_ratio]
            
            ratio_df = pd.DataFrame(ratio_list, columns=['Gene','Guide', 'Ratio'])
            
            os.chdir(csv_folder)
            ratio_df.to_csv(csv_folder+'/ratios/' +gene_name+ '_CelFi_ratios.csv',index=False)


        elif(num_of_guides == 1):
            firstguide_first_time = int(input("\nEnter row number of first time point: "))
            firstguide_last_time = int(input("\nEnter row number of last time point: "))
            
            #firstguide_first_time = 2
            #firstguide_last_time = 3
            
            guide_name.append(str(oof_df.at[firstguide_first_time,'Sample']).split(' ')[0])
            
            ratio.append(round(oof_df.at[firstguide_last_time, 'Out-of-frame']/oof_df.at[firstguide_first_time, 'Out-of-frame'],2))
            
            print(f"Guide: {guide_name[0]}")
            guide_name_good = input("Please confirm guide names, y or n: ").upper()
            
            if guide_name_good == 'N' or guide_name_good == 'NO':
                guide_name.clear()
                guide_name.append(input("Guide name: "))

            print()
            print(f"{guide_name[0]} fitness ratio: {ratio[0]}")

            
            fig, ax = plt.subplots()
            #makes the individual guide bar plot a separate entity for easier editing
            
            x_pos = [0,1,2]
            
            bars = ax.bar(x_pos,
                        [0,ratio[0],0], 
                        color='#16A085', 
                        edgecolor='#17202A',
                        width=0.5,
                        )
            #lazy way of compiling ratios for csv export
            g1_ratio = [gene_name, guide_name[0], ratio[0]]

            
            ratio_list = [g1_ratio]
            
            ratio_df = pd.DataFrame(ratio_list, columns=['Gene','Guide', 'Ratio'])
            
            os.chdir(csv_folder)
            ratio_df.to_csv(csv_folder+'/ratios/' +gene_name+ '_CelFi_ratios.csv',index=False)                  
            label_bar(bars)
            plt.xticks(range(1,2),guide_name)
        
        plt.ylim(ymax=1.2, ymin=0)
        plt.ylabel('Fitness Ratio')
        plt.title(chart_title ,y=1.05)
        ax.spines[['top','right']].set_visible(False)

        plt.show()
    
    
def _compile_ratios():
    
    ratio_dir = working_dir+"/_cell_fitness_folder/csv_files/ratios"
    os.chdir(ratio_dir)
    ratio_list=[]
    compiled_ratios = pd.DataFrame()
    for csv in os.scandir(ratio_dir):
        ratio_list.append(csv)

    for ratio in ratio_list:
        tmp = pd.read_csv(ratio)
        compiled_ratios = compiled_ratios.append(tmp, ignore_index=True)
        
    #compiled_ratios.drop(columns=['Unnamed'])
    print(compiled_ratios.head())
    
    print(os.getcwd())
    input()
    compiled_ratios.to_csv('all_ratios.csv',index=False)

#_run_all_indels()


working_dir, cagenum_dirs,cage_num = _input(working_dir, NGS_DATE)
        
target_csv,graph_dir = _rename_csv(working_dir, cagenum_dirs)


_make_plot()



#_compile_ratios()




print()
print("Program completed.")

