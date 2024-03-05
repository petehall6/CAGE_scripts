from Bio import SeqIO

input_file = "ENTER YOUR INPUTFILE NAME HERE. name something other than gene.gb"

input_handle = open(input_file,'r')
output_handle = open("gene.gbk","w+")

seqs = list(SeqIO.parse(input_handle, "genbank"))

write_out = SeqIO.write(seqs, output_handle, "genbank")

output_handle.close()
input_handle.close()