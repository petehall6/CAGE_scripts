
from ssODN_for_deletions import Main
from loxp_model import lox_p_Main
from loxP_ssODN_generator_v3 import create_loxP_from_core_projects_dir
import sys

sys.path.insert(1, "/research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/python/homology_arm_primer_finder")
from Junction_Primer_Generator_400bp_HA import jxn_primer_gen

def choose_program():
    print("\n\n                 ********* Donor Design Programs: *********\n")
    print("\n1: Junction_Primer_Generator (designs primers surrounding left and right homology arms)\
        \n2: ssODN_for_deletions (designs bridging ssODNs for deletions)\
        \n3: loxP_ssODN_generator (designs loxP ssODNs)\
        \n4: loxP_model (combines 2 CAGE files, creating a gbk with all gRNAs, a gbk with a deletion, and (LOXP IS NOT CORRECT")

    choice = str(input("\nWhich program would you like to run? (type 1, 2, or 3 and hit enter): "))

    
    #junction primer
    if choice == "1":
        print("\n  * Running Junction_Primer_Generator...\n\n")
        jxn_primer_gen()
    #ssODN deletion generator
    if choice == "2":
        print("\n  * Running ssODN_for_deletions...\n\n")
        # TODO: RENAME
        Main()
    #loxP ssODN generator
    if choice == "3":
        print("\n  * Running loxP_ssODN_generator...\n\n")

        
        create_loxP_from_core_projects_dir()
    if choice == "4":
        print("\n  * Running loxP_model...\n\n")
        lox_p_Main()

def main():
    loop = "Y"
    while loop == "Y":
        choose_program()
        loop = input("Would you like to run another design program? (y/n): ").upper()
    return

if __name__ == '__main__':
    main()
