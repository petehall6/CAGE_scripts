#7/6/2023 PMH

#This script plots Fitness Assay NGS Data

import os
import pandas as pd
import shutil
import win32com.client
import matplotlib as mpl
import matplotlib.pyplot as plt




NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"

def _input():
    #get input
    #run_date = "062023"
    run_date = input("Enter run date: ").strip()
    working_dir = os.path.join(NGS_DIR, run_date, "joined").replace("\\","/")

    #check if run_date is correct/folder exsists
    while os.path.isdir(working_dir) == True:
        print("Success!  Run date found.")
        break
        
    else:
        print("Run date not found.  Please reenter run date.")
        run_date = input("Enter run date: ").strip()
        working_dir = os.path.join(NGS_DIR, run_date, "joined").replace("\\","/")
        
    
    cagenum_dirs =[]
    #check project number/cage number.  If there anything in the cagenum/projectnum dirs list loop breaks and contiues to emailer
    while len(cagenum_dirs) == 0:
        cage_num = input("Please enter project number: ").upper().strip()
        #cage_num = "CAGE9999"
        print("Searching for project folder.  This may take a moment.")
        cagenum_dirs = [dir.name for dir in os.scandir(working_dir) if dir.is_dir() if cage_num in dir.name]
        print(f"Number of projects found: {len(cagenum_dirs)}")
        #empty list evaluate to false
        if not cagenum_dirs:
            print(f"Did not find project number {cage_num}.")
                    
    return working_dir,cagenum_dirs,cage_num

def _rename_csv(working_dir, cagenum_dirs):
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
        del raw_df
        
        print("\nNaming complete\n\n")

    print(target_csv)
    return(target_csv,graph_dir)

def _make_plot(target_csv, graph_dir):
    
    def label_bar(bars):
            for bar in bars:
                height = bar.get_height()
                if height > 0:
                    ax.text(bar.get_x() + bar.get_width()/2.0, 1.0*height,
                            f'{height}',
                            ha='center', va='bottom')
    
    
    guide_names = []
    guide_scores = []
    
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
        #guide_names = ['Guide One','Guide Two']
        #firstguide_first_time = 2
        #firstguide_last_time = 3
        #secondguide_first_time = 4
        #secondguide_last_time = 5
        
        
        firstguide_first_time = int(input("\nEnter row number of first time point for first guide: "))
        firstguide_last_time = int(input("\nEnter row number of last time point for first guide: "))
        
        secondguide_first_time = int(input("\nEnter row number of first time point for second guide: "))
        secondguide_last_time = int(input("\nEnter row number of last time point for second guide: "))
        
        
        guide_names.append(str(oof_df.at[firstguide_first_time,'Sample']).split('.')[1].split(' ')[0])
        guide_names.append(str(oof_df.at[secondguide_first_time,'Sample']).split('.')[1].split(' ')[0])
        
        
        print(f"1st guide: {guide_names[0]}")
        print(f"2nd guide: {guide_names[1]}")
        
        guide_names_good = input("Please confirm guide names, y or n: ").upper()
        
        if guide_names_good == 'N' or guide_names_good == 'NO':
            guide_names.clear()
            guide_names.append(input("1st guide name: "))
            guide_names.append(input("2nd guide name: "))
        
        guide_scores.append(round(oof_df.at[firstguide_last_time, 'Out-of-frame']/oof_df.at[firstguide_first_time, 'Out-of-frame'],2))
        
        guide_scores.append(round(oof_df.at[secondguide_last_time, 'Out-of-frame']/oof_df.at[secondguide_first_time, 'Out-of-frame'],2))
        
            
        print()
        print()
        print(f"{guide_names[0]} fitness ratio: {guide_scores[0]}")
        print()
        print()
        print(f"{guide_names[1]} fitness ratio: {guide_scores[1]}")
        
        #sets overall graph as subplot for easier editing
        fig, ax = plt.subplots()
        #makes the individual guide bar plot a separate entity for easier editing
        
        
        
        #Okay, so this is a super ugly hacky way of doing this....
        #x_pos sets the column number for each bar plot.  In the case of two guides there are two empty columns flanking the real guides
        x_pos=[0,1,2,3]
        bars = ax.bar(x_pos,
                      # sets values for empty columns as 0 and explicitly calls guide score indices
                      [0,guide_scores[0],guide_scores[1],0], 
                      color='#16A085', 
                      edgecolor='#17202A',
                      width=0.5,
                      )


        label_bar(bars)
        #uses a short range of 1,2 to get two middle x_pos columns and uses guide_names list to fill in names
        plt.xticks(range(1,3), guide_names)
        
        #lazy way of compiling ratios for csv export
        g1_ratio = [gene_name, guide_names[0], guide_scores[0]]
        g2_ratio = [gene_name, guide_names[1], guide_scores[1]]
        
        ratio_list = [g1_ratio,g2_ratio]
        
        ratio_df = pd.DataFrame(ratio_list, columns=['Gene','Guide', 'Ratio'])
        
        os.chdir(graph_dir)
        ratio_df.to_csv(gene_name+ '_CelFi_ratios.csv')


    elif(num_of_guides == 1):
        firstguide_first_time = int(input("\nEnter row number of first time point: "))
        firstguide_last_time = int(input("\nEnter row number of last time point: "))
        
        #firstguide_first_time = 2
        #firstguide_last_time = 3
        
        guide_names.append(str(oof_df.at[firstguide_first_time,'Sample']).split('.')[1].split(' ')[0])
        
        guide_scores.append(round(oof_df.at[firstguide_last_time, 'Out-of-frame']/oof_df.at[firstguide_first_time, 'Out-of-frame'],2))
        
        print(f"Guide: {guide_names[0]}")
        guide_names_good = input("Please confirm guide names, y or n: ").upper()
        
        if guide_names_good == 'N' or guide_names_good == 'NO':
            guide_names.clear()
            guide_names.append(input("Guide name: "))

        print()
        print(f"{guide_names[0]} fitness ratio: {guide_scores[0]}")
        
        
        
        
        
        fig, ax = plt.subplots()
        #makes the individual guide bar plot a separate entity for easier editing
        
        x_pos = [0,1,2]
        
        bars = ax.bar(x_pos,
                      [0,guide_scores[0],0], 
                      color='#16A085', 
                      edgecolor='#17202A',
                      width=0.5,
                      )
        #lazy way of compiling ratios for csv export
        g1_ratio = [gene_name, guide_names[0], guide_scores[0]]

        
        ratio_list = [g1_ratio]
        
        ratio_df = pd.DataFrame(ratio_list, columns=['Gene','Guide', 'Ratio'])
        
        os.chdir(graph_dir)
        ratio_df.to_csv(gene_name+ '_CelFi_ratios.csv')                    
        label_bar(bars)
        plt.xticks(range(1,2),guide_names)
    plt.ylabel('Fitness Ratio')
    plt.title(chart_title ,y=1.05)
    ax.spines[['top','right']].set_visible(False)

    plt.show()
    
    




    
    

working_dir, cagenum_dirs,cage_num = _input()
        
target_csv,graph_dir = _rename_csv(working_dir, cagenum_dirs)


#target_csv = r"Z:/ResearchHome/Groups/millergrp/home/common/NGS/062023/joined\\CAGE9999TEST\\CAGE9999_hPSMD12_F2_R2_all_indels.csv"
#graph_dir = r"Z:/ResearchHome/Groups/millergrp/home/common/NGS/062023/joined\\CAGE9999TEST\\graphs"

_make_plot(target_csv, graph_dir)


print()
print("Program completed.")

