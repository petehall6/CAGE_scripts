import os
import glob
import time
#way to scroll through NGS run date, CAGE#, well# and pull fastq files for automated zipping and emailing.
NGS_DIR = "Z:/ResearchHome/Groups/millergrp/home/common/NGS"

os.chdir(NGS_DIR)
print(os.getcwd())

run_date = '111423'

working_dir = os.chdir(os.path.join(NGS_DIR,run_date,"joined"))
print(os.getcwd()) 

#search NGS folder and find the results csv
cage_nums = ['CAGE104']
cage_folders = []
start_time = time.time()
for num in cage_nums:
    print(num)
    #find all matches with CAGe# and only return directories
    cage_folders.append([name for name in glob.glob(f"{num}*") if os.path.isdir(name)])
end_time = time.time()
print(cage_folders)
print(f"Run time: {end_time-start_time}")

#