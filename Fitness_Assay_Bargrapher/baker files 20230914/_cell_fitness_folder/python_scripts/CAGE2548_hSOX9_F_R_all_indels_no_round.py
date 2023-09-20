import os,sys,inspect
import argparse

def main():
    ID = str(os.path.basename(__file__)).split('.py')[0]
    the_Seq = str.upper('CAGCCCCACCATGTCCGAGGACTCCGCGGGCTCGCCCTGCCCGTCGGGCTCCGGCTCGGACACCGAGAACACGCGGCCCCAGGAGAACACGTTCCCCAAGGGCGAGCCCGATCTGAAGAAGGAGAGCGAGGAGGACAAGTTCCCCGTGTGCATCCGCGAGGCGGTCAGCCAGGTGCTCAAAGGCTACGACTGGACGCTGGTGCCCATGCCGGTGCGCGTCAACGGCTCCAGCAAGAACAAGCCGCACGTCAAGCGGCCCATGAACGCCTTCATGGTGTGGGCGCAGGCGGCGCGCAGGAAGCTCGCGGACCAGTACCCGCACTTGCACAACGCCGAGCTCAGCAAGACGCTGGGCAAGCTCTGGAGGTAGGACCCGGC')
    the_Seq_start = str.upper('TGCCCGTCGGGCT')
    the_Seq_end = str.upper('AACGCCGAGCTCA')
    fastq_files = 'Miller-Plate[1][3][1-5]*.fastq'
    test_list = [str("g19"), str.upper("GGCGGCGCGCAGGAAGCTCGCGG"), str("g36"), str.upper("CCCAAGGGCGAGCCCGATCTGAA")
				]



    search_fastq(ID,the_Seq,the_Seq_start,the_Seq_end,fastq_files,test_list)





if __name__ =='__main__':		# use this if you want to include modules from a subfolder
    cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(os.path.split(inspect.getfile( inspect.currentframe() ))[0], r'Z:\ResearchHome\Groups\millergrp\home\common\Python\CAGE_Programs\NGS_one_main_program')))
    if cmd_subfolder not in sys.path:
        sys.path.insert(0, cmd_subfolder)
        from clone_check_main_all_indels_no_round_v3 import search_fastq
    main()
