import os
import glob
ngs_dir = r'Z:\ResearchHome\Groups\millergrp\home\common\NGS'

out_file = 'all_indels_list.txt'

dir_list = [ '071624',
            '071124',
            '070924',
            '070224',
            '062724',
            '062524',
            '061824',
            '061324']

all_indels_list=[]
print("scanning")
for dir in dir_list:
    target = os.path.join(ngs_dir,dir,'joined')
    for obj in os.scandir(target):
        if obj.is_dir() and obj.name.endswith("all_indels"):
            all_indels_list.append([dir,obj.name])
    



print(all_indels_list)
print('saving')

with open(out_file, 'w+') as f:
    for indel in all_indels_list:
        f.write(f"{indel}\n")
    f.close()
    

