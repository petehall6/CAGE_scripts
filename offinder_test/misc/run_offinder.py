import os
import subprocess
import re

#* Change input name if needed
#* should be formatted as: sequence plus PAM 'space' mismatch number 'space' guide name.  ie:  AGGGAGGAAGGGAGCCTCAANGG 3 GUIDE.gNumber
input_guides = 'input_CasOffinder.txt'

print("loading cas-offinder/2.4.1.  Ignore the note below about it only works with in the gpu queue")
os.system('module load cas-offinder')

offinder_call = f'bsub -P Cas_offinder -q rhel8_gpu_short -gpu "num=1/host" -R a100 -M 8000 -o out_CasOffinder.out -e err_CasOffinder.err -J cas_offinder "hostname & module load cas-offinder && cas-offinder {input_guides} G output_CasOffinder.txt"'
                


Long_sout = subprocess.run([offinder_call], shell=True, stdout=subprocess.PIPE).stdout.decode('utf-8')    
jid_Long9 = re.split("<|>", Long_sout)[1]
if jid_Long9.isnumeric():
    print(f"JobID-{jid_Long9} for cas_offinder has been successfully submitted. Please wait.")
    print("System will notify when the job is complete")
else:
    print("Failed in submitting job for cas-offinder OTA for Cas9Long")
os.system(f"bwait -w 'ended(cas_offinder)'")
print("Jobs done")