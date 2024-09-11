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
python library_from_crispick.py $1 $2 $3 $4 $5 $6 5000