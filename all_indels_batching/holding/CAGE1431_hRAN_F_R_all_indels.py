
#TODO make multi_range template


import os,sys,inspect
import argparse

def main():
    ID = str(os.path.basename(__file__)).split('.py')[0]
    the_Seq = str.upper('aggtgactcgttgaggtctaagactataactagttgaatcctagtctctggttcttagcattcttcacttgtgcagtaaactgagtgtactaattcccacaaatgtttcttcagcccagtgtgccatcataatgtttgatgtaacatcgagagttacttacaagaatgtgcctaactggcatagagatctggtacgagtgtgtgaaaacatccccattgtgttgtgtggcaacaaagtggatattaaggacaggaaagtgaaggcgaaatccattgtcttccaccg')
    the_Seq_start = str.upper('agtctctggttct')
    the_Seq_end = str.upper('gacaggaaagt')
    #fastq_files = *.fastq
    test_list = [str("g3"), str.upper("aacatccccattgtgttgtgtgg"), str("g10"), str.upper("cccattgtgttgtgtggcaacaa")]
    
    
    glob_range = ['Miller-Plate25*.fastq', 'Miller-Plate26*.fastq', 'Miller-Plate30*.fastq']
    
    search_fastq(ID,the_Seq,the_Seq_start,the_Seq_end,glob_range,test_list)
    
    
if __name__ =='__main__':		# use this if you want to include modules from a subfolder
        cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile(inspect.currentframe() ))[0], r'Z:\ResearchHome\Groups\millergrp\home\common\Python\CAGE_Programs\NGS_one_main_program')))
        if cmd_subfolder not in sys.path:
            sys.path.insert(0, cmd_subfolder)
            from clone_check_main_multi_all_indels import search_fastq
        main()
    