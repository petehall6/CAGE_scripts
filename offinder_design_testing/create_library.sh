#!/bash/bin

module load bowtie
module load primer3
#for new rhel8 servers
module load cas-offinder/2.4.1
#module list
module list

current_dir=$(pwd)

python accessory_files/library_design.py  --cwd "$current_dir" "$@" 