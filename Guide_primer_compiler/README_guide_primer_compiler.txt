README_guide_primer_compiler

This program will parse through all of the core project folders to find all projects that have the same gene.  It will then create a blank genbank file and add all the primers, guides and CDS features from all of the associated projects.

To Use:

1) Enter gene name
2) Select species by entering h for human or m for mouse

NOTE:

If you get a message similar to the following:

	Failed to add feature from project @ : Z:\ResearchHome\Groups\millergrp\home\common\CORE PROJECTS\hWHAMM-Schultz-Cherry-CAGE2265\WHAMM_mod_NGS.gbk)
	No note added to feature: OrderedDict([('blurb', ['CAGE2265.WHAMM.Cas9.g3'])])


This can happen.  The program is looking for misc_features that have a "note" inside of the genbank file.  The error message will direct you to the project and the feature.  You can simply add the correct note to add the feature or ignore it if the feature is not needed, ie a label feature was added to the genbank file but not necssary for sequence information.  

If you want to add the feature that was not added:
	1) Open the genbank file
	2) Navigate to features and find the problem feature
	3) Double click the name
	3) Enter information in the bottom text box and make sure the the drop down box is "/note"
	4) Rerun script