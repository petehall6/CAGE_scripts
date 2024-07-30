#!/bash/bin

# Tests designs file in current conda environment
# NOTE: DO NOT 'module load' PYTHON UP/DOWNSTREAM.
#      THAT WILL PREVENT THIS RUNNING IN CONDA ENV.
module load bowtie
module load primer3
#for new rhel8 servers
module load cas-offinder/2.4.1-rhel8

module list

echo Compiling OTA scores and generating tables
#module list
python /research_jude/rgs01_jude/groups/millergrp/projects/CrisprDesigns/common/screens/programs/offinder_test/from_a_guide_list/get_ota_from_guide_list.py 5000

