#!/usr/bin/env python

import os


def check_projects(gene_name, species):
    #    print 'Checking to see if we have done any similar projects to {}, {}'.format(species,gene_name)

    if species == "mouse":
        gene_name = "m" + gene_name
    elif species == "human":
        gene_name = "h" + gene_name
    else:
        pass

    core_projects_dir = (
        "/research_jude/rgs01_jude/groups/millergrp/home/common/CORE PROJECTS"
    )
    in_progress_dir = os.path.join(core_projects_dir, "Projects_in_progress_ppts")
    folders_to_check = [core_projects_dir, in_progress_dir]

    ppt_list = []

    def _check_subdirs_for_gene(path):
        for f in os.scandir(path):
            if f.is_dir():
                if gene_name.upper() in f.name.upper():
                    ppt_list.append(str(f.path))

    for parent_folder_path in folders_to_check:
        _check_subdirs_for_gene(parent_folder_path)

    if len(ppt_list) == 0:
        pass
        # print 'No projects with that gene'
    else:
        print(f"\n\n btw, did you know we have done projects involving {gene_name} in a {species}?")
        print("These files may be useful:")
        for x in ppt_list:
            print(f"   * {x}")


def main():
    # Testing defaults
    check_projects("HBB", "human")


if __name__ == "__main__":
    main()