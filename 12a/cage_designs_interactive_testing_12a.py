#!/usr/bin/env python
"""
This script will break if it is not in the same directory as the cage_designs.py script
"""

import cmd

# import logging
import os

# import pdb
# import subprocess
# import sys
import pandas as pd
import numpy as np
from datetime import datetime

# import datetime
from getpass import getuser
import KO_blacklist_check
import srmFetcher_objective
import crisprSummarySRMConverter
from cage_designs_testing_12a import (
    BIN_DIR,
    DESIGNS_DIR,
    SUPPORTED_SPECIES,
    VECTOR_CHOICES,
    Z_DESIGNS_DIR,
    Z_DESIGNS_FOR_CLUSTER,
    EntrezGeneClient,
    GuideRNA,
    make_project_dir,
    SRMJoiner
)

import ess_check
import CopyNum


# bash_var = 5000


def format_gene_name(gene_name, species):
    if species.lower() == "human":
        gene_name = str(gene_name).upper()
    else:
        gene_name = str(gene_name).title()
    return gene_name


class CrisprDesign(cmd.Cmd):
    menu = """
                    MAIN MENU
    +-------------+-----------+-----------+
    |   DESIGNS   |   TOOLS   |  OPTIONS  |
    +-------------+-----------+-----------+
    |             |           |           |
    |    init*    |  joiner   |  size     |
    |     OR      |           |           |
    |   manual    |  repick   |  help     |
    |             |           |           |
    |    grna     |  convert  |  menu     |
    |             |           |           |
    |   primer    |  copies   | allfeats  |
    |     OR      |           |           |
    |   pfuzzy    |           |           |
    |             |           |           |
    |             |           |           |
    +-------------+-----------+-----------+
    * "init" will only show designs assigned to you.
        "init all" will show current designs in SRM.

            Not sure what a command does?
        Type "help [command]" for more info!
"""
    intro = (
        f"{menu}"
        + "\nPress Ctrl-D to exit a command.\nPress Ctrl-C to exit program"
        + "\nTo run shell commands type ! then command."
        + "\nTo view most recent files type !ls -lt\n"
        # + "\n"
    )
    prompt = """
    To display the main menu again, run the 'menu' command!
    
    Enter a command: """

    def __init__(self):
        self.BIN_DIR = os.path.dirname(os.path.abspath(__file__))
        # will be set by init command or passed by user in subsequent commands
        self.PROJECT_DIR = None
        self.gle = 5000
        self.allfeats = False
        cmd.Cmd.__init__(self)

    def do_menu(self, line):
        """
    * MENU
        Print the main menu again.
        """
        try:
            print(self.menu)
        except EOFError:
            return None

    def do_joiner(self, line):
        """
    * JOINER
        Run SRM2 multi-txt file joiner
        """
        try:
            print()
            SRMJoiner()
        except EOFError:
            return None

        # system_call = os.path.join(self.BIN_DIR, "cage_designs.py srmjoiner")
        # print(system_call)
        # try:
        #     os.system(system_call)
        #     return None
        # except EOFError:
        #     return None

    def do_repick(self, line):
        """
    * REPICK
        Allows you to re-pick gRNAs from existing CRISPR_Summary file.
        """
        crisprSummarySRMConverter.repick(
            os.path.join(DESIGNS_DIR, self.prompt_project_dir())
        )
        return None

    def do_convert(self, line):
        """
    * CONVERT
        Run SRM2 Converter for a single TXT file.
        Helpful if you had to change picks and need to re-upload.
        """
        try:
            print()
            return crisprSummarySRMConverter.singleConverter()
        except EOFError:
            return None

    def do_allfeats(self, line):
        print("Would you like to add all gRNA features to gbk?")
        try:
            response = input("  --> ").lower().strip()
            if response in ["y", "yes"]:
                print("\nConsider it done!\n")
                self.allfeats = True
            elif response in ["n", "no"]:
                print("\nConsider it done!\n")
                self.allfeats = False
            else:
                print(
                    "I don't know what your response means, so I'm returning to the main menu."
                )
            return None
        except EOFError:
            return None

    def do_size(self, line):
        """
    * SIZE
        Sets the size (kilobase) of the extension on either side of the sequence.
        """
        print(f"\nCurrent size is: {self.gle}")
        print("Would you like to change it? (y/n)\n")
        try:
            response = input("  --> ").lower().strip()
            if response in ["y", "yes"]:
                print("\nEnter desired number of kb. (5, 500, etc)")
                kb_choice = input("  --> ").strip()
                try:
                    kb_choice = int(kb_choice)
                    self.gle = kb_choice * 1000
                except:
                    print(
                        "Could not convert kb choice to an integer. Please set this option again before continuing."
                    )
            return None
        except EOFError:
            return None

    def do_manual(self, line):
        """
    * MANUAL
        Manually enter information to initialize a design.
        Helpful for projects not in SRM2!
        """
        GENBANK_LENGTH_EXTENSION = self.gle
        args = ["gene", "customer", "experiment"]
        system_call = ""
        try:
            values = self.prompt_user_for_cmdline_args(args)
            values["species"] = self.prompt_user_enforce(
                "species",
                "Choose from [ %s ]" % ", ".join(SUPPORTED_SPECIES),
                lambda x: x in SUPPORTED_SPECIES,
            )
            # values['vector'] = self.prompt_user_enforce(
            #     'vector', 'Choose from [ %s ]. Type help to see description' % ', '.join(VECTOR_CHOICES),
            #     lambda x: x in VECTOR_CHOICES)
        except EOFError:
            return None

        values["GENBANK_LENGTH_EXTENSION"] = self.gle
        values_copy = values.copy()
        entrezClient = EntrezGeneClient()
        values_copy["gene"] = entrezClient.gene_token_to_name(values["gene"])
        values_copy["gene"] = format_gene_name(
            values_copy["gene"], values_copy["species"]
        )
        values["gene"] = values_copy["gene"]

        system_call = os.path.join(self.BIN_DIR, "cage_designs_testing_12a.py")
        system_call += " init -g {gene} -c {customer} -e {experiment} -sp {species} -gle {GENBANK_LENGTH_EXTENSION}"
        system_call = system_call.format(**values)
        self.PROJECT_DIR = make_project_dir(**values_copy)
        os.system(system_call)
        return None

    def do_init(self, line):
        """
    * INIT
        Setup a project directory inside the current working directory
        Defaults to "selfish" fetch from SRM2, showing only designs assigned to you.
        To show all SRM2 designs, run 'init all'
        To initialize a design manually, run 'manual'
        """
        args = ["gene", "customer", "experiment"]
        system_call = ""
        GENBANK_LENGTH_EXTENSION = self.gle

        def chooseExperiment(
            gene,
            customer,
            experiment,
            species,
            listPath=os.path.join(Z_DESIGNS_FOR_CLUSTER, "cage_projects_list.xlsx"),
        ):
            # Latest experiment is CAGE###.
            # Press y to use this number, or type custom experiment number now.
            def _read_current_xlsx():
                """
                Reads Excel file listing all CAGE projects in as a DataFrame.
                Updates to include integers for the project number if missing.
                Returns resultant DataFrame.
                """

                def __create_cage_num(series):
                    series.cage_number = series.expmt.replace("CAGE", "")
                    return series

                expmt_df = pd.read_excel(
                    listPath, sheet_name="projects", parse_dates=["date"]
                )
                cage_df = expmt_df[expmt_df["expmt"].str.contains("CAGE") == True]
                nan_df = cage_df[np.isnan(cage_df["cage_number"])]
                nan_df = nan_df.sort_values(["cage_number"])
                nan_df = nan_df.apply(
                    __create_cage_num, axis=1, result_type=None
                ).astype({"cage_number": "int64"}, errors="ignore")
                expmt_df.update(nan_df["cage_number"])
                return expmt_df

            def _append_and_save(target_df, experiment_title, cage_number_value):
                target_df = target_df.append(
                    {
                        "expmt": experiment_title,
                        "cage_number": cage_number_value,
                        "date": datetime.now(),
                        "designer": getuser(),
                        "gene": gene,
                        "customer": customer,
                        "experiment": experiment,
                        "species": species,
                        "notes": "",
                    },
                    ignore_index=True,
                )
                target_df.to_excel(
                    listPath, sheet_name="projects", header=True, index=False,
                )
                return target_df

            current_df = _read_current_xlsx()
            new_number = int(current_df["cage_number"].max()) + 1
            print(f'Latest CAGE Project number would be {f"CAGE{new_number}"}.')
            print("Enter y to proceed, or enter custom project number now.")
            try:
                response = input("  --> ").strip().upper()
                if response == "Y":
                    results = f"CAGE{new_number}"
                    current_df = _append_and_save(current_df, results, new_number)
                    return results
                else:
                    results = response.strip()
                    current_df = _append_and_save(current_df, results, "")
                    return results
            except EOFError:
                return None

        selfish = True
        if "all" in line:
            selfish = False
        try:
            (gene, customer, experiment, species, cell_line, objective) = srmFetcher_objective.main(selfish)
        except TypeError:
            return None
        except Exception as error:
            print(
                "\nFailed in design fetcher! Please take a screenshot and notify programming team."
            )
            print(error)
            return None
        if cell_line:
            KO_blacklist_check.blacklist_check(cell_line)
        if gene and cell_line:
            ess_check.cellLineEssential_Dependency(gene, cell_line)
            
        
        CopyNum.getCopyNumberInLine(cell_line, gene, objective)
        gene = gene.strip()
        gene = format_gene_name(gene, species)
        
        system_call = os.path.join(self.BIN_DIR, "cage_designs_testing_12a.py")
        system_call += f" init -g {gene} -c {customer} -e {experiment} -sp {species} -gle {GENBANK_LENGTH_EXTENSION}"
        projectDir_kwargs = {
            "gene": gene,
            "customer": customer,
            "experiment": experiment,
            "species": species,
        }
        self.PROJECT_DIR = make_project_dir(**projectDir_kwargs)

        # print()
        # print("System call is: ")
        # print(f"     {system_call}")
        # print()
        os.system(system_call)
        return None

    def do_grna(self, line):
        """
    * GRNA
        To pick gRNAs in a project directory that has already had gRNAs picked, type "grna force"
        """
        args = ["sequence"]

        try:
            values = {}
            values["project_directory"] = os.path.basename(self.prompt_project_dir())
            print("\nVECTOR CHOICES:")
            print("  SP78 (use instead of SNP23!): for Mammalian U6")
            print("  SYN: Synthetic gRNA")
            print("  SNP23 (aka SM115): for Mammalian U6")
            print("  SNP21 (aka SM102): for Zebrafish T7")
            print("  SNP22 (aka SM104): for pX330 single vector gRNA + Cas9")
            print("  SM387: Lenti-CRISPR")
            print("  SNP36 (aka SP82):  Lenti-CRISPRv2")
            print("  SNP329:  Lenti-CRISPRv2 + mCherry")
            print("  SNP112: (Lenti-CRISPRv2 w eGFP)")
            print(
                "  SNP27 (aka SM388): pX458 single vector gRNA +Cas9-T2A-eGFP (AppliedStemCell)"
            )
            print("  SNP38 (aka SS32, PC323) pX459 single vector gRNA +Cas9-T2A-Puro")
            print(
                "  SNP28 (aka SM438): lentiGuide-Puro  single vector gRNA + EF1a-Puro (Milbrandt)"
            )
            print("  SP16: lentiGuide-eGFP, single vector gRNA that expresses eGFP")
            print("  SP76 : Sanger vector")
            print("  PC291 : Activation vector with Puro selection")
            print()
            values["vector"] = self.prompt_user_enforce(
                "vector",
                "",
                # "Choose from [ %s ]" % ", ".join(VECTOR_CHOICES),
                lambda x: x.upper() in VECTOR_CHOICES,
            )
            values.update(self.prompt_user_for_cmdline_args(args))

        except EOFError:
            return None

        system_call = (
            os.path.join(self.BIN_DIR, "cage_designs_testing_12a.py")
            + ' grna -p {project_directory} -v {vector} -s "{sequence}"'
        )
        system_call += " --force" if "force" in line else ""
        system_call = system_call.format(**values)

        print(f"\nSystem Call:\n    {system_call}")
        os.system(system_call.format(**values))
        # subprocess.check_call( system_call.format(**values).split())

    def do_primer(self, line):
        """
    * PRIMER
        Design primers for desired gRNAs.
        """
        args = ["guide_rna_IDs"]

        try:
            values = {}
            values["project_directory"] = os.path.basename(self.prompt_project_dir())
            # values["fuzzy"] = self.prompt_user_enforce(
            #     "Fuzzy match control string ie {s<=1}",
            #     "Default = None",
            #     lambda x: True if x is None else GuideRNA.validate_fuzzy_match(x),
            # )

            values.update(self.prompt_user_for_cmdline_args(args))

        except EOFError:
            return None

        system_call = (
            self.BIN_DIR
            + os.path.sep
            #TODO changed file
            + "cage_designs_testing_12a.py primer -p {project_directory} "
            + "-i {guide_rna_IDs}"
        )
        # system_call += ' --fuzzy "{fuzzy}"' if values["fuzzy"] else ""
        system_call = system_call.format(**values)
        # print(system_call)
        os.system(system_call)

        # TODO: Make spit out link to "View Details" on design page
        # LINK_TEMPLATE =
        # https://srm.stjude.org/srm/lims/taskDataEntry.html?serviceId=120&navigationId=
        # Worklog{2817263}
        # Order{206269}
        # OrderEavObject{6225753}
        # Subject{2072323}
        # SubjectEavObject{6225752}
        # Task{3158714}
        # TaskEavObject{6225806}
        # &businessActivityId=32033

    def do_pfuzzy(self, line):
        """
    * PFUZZY
        Design primers using "fuzzy matching."
        This allows you to design primers for the regular gene
            even when your target sequence was mutant.
        
        SYNTAX:
        * Curly braces are required
        * Followed by one of following: i, d, s
        ** (Short for insertion, deletion, substitution.)
        * Followed by logical operator (<, >, =, <=, ...)
        * And your desired threshhold.

        i.e. {s<=1} would allow for up to one substitution.

        NOTE: Multiple fuzzy constraints can be used at once.
            Separate multiple by a comma (,) i.e. {i=0,s<=1}
        """
        args = ["guide_rna_IDs"]

        try:
            values = {}
            values["project_directory"] = os.path.basename(self.prompt_project_dir())
            values["fuzzy"] = self.prompt_user_enforce(
                "Fuzzy match control string ie {s<=1}",
                "Default = None",
                lambda x: True if x is None else GuideRNA.validate_fuzzy_match(x),
            )

            values.update(self.prompt_user_for_cmdline_args(args))

        except EOFError:
            return None

        system_call = (
            self.BIN_DIR
            + os.path.sep
            + "cage_designs_testing_12a.py primer -p {project_directory} "
            + "-i {guide_rna_IDs}"
        )
        system_call += ' --fuzzy "{fuzzy}"' if values["fuzzy"] else ""
        system_call = system_call.format(**values)
        # print(system_call)
        os.system(system_call)
        return None

    def do_copies(self, line):
        """
    * COPIES
        To check copy number of a gene.  Stand-alone format"
        """
        try:
            print("\n\nCheck copy number of gene \n\n")

            CopyNum.getCopyNumberTools()

        except EOFError:
            return None
    def do_shell(self, line):
        """
    * SHELL
        Run a shell command
        """
        print("running shell command:", line)
        output = os.popen(line).read()
        print(output)

    def emptyline(self):
        pass

    def prompt_project_dir(self):
        if self.PROJECT_DIR == None:
            self.PROJECT_DIR = DESIGNS_DIR
            return self.prompt_project_dir()
        else:
            project_dir = self.prompt_user_enforce(
                "Project Directory",
                "Default = %s" % os.path.basename(self.PROJECT_DIR),
                lambda x: os.path.exists(os.path.join(DESIGNS_DIR, x)),
                default=self.PROJECT_DIR,
            )

            self.PROJECT_DIR = os.path.join(DESIGNS_DIR, project_dir)
            return project_dir

    def prompt_user_for_cmdline_args(self, args):
        values = {}
        try:
            for arg in args:
                choice = str(input("Enter %s -> " % arg))
                values[arg] = choice.strip()

            return values
        except EOFError:
            return None

    def prompt_user_enforce(self, arg, prompt_append, condition, default=None):
        """
        arg is name of variable, prompt_append is string to add to end of prompt
        condition is a function that must be true before returning
        default is value used for an empty line
        """
        try:
            while True:
                user_input = input("Enter %s : %s -> " % (arg, prompt_append)).strip()

                value = user_input if user_input else default

                if condition(value):
                    return value
                else:
                    print("Input %s not valid for %s " % (value, arg))
        except EOFError:
            return None


if __name__ == "__main__":
    try:
        CrisprDesign().cmdloop()
    except KeyboardInterrupt:
        print("\n\n")
        pass
