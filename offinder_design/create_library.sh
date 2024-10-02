#!/bash/bin

module load bowtie
module load primer3
#for new rhel8 servers
module load cas-offinder/2.4.1-rhel8

module list

#module list
python library_design.py "$@"
