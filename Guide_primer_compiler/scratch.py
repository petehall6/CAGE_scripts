from Bio import SeqIO, SeqFeature, SeqRecord, GenBank
import os


template_gbk = "gene.gbk"
found_features = []

output_handle = open("new_gbk.gbk", "w")

for record in SeqIO.parse("ADA2_mod_NGS.gbk", "genbank"):
    for feature in record.features:
        if feature.type == "primer_bind":
            found_features.append(feature)
            print(feature)
        
        if feature.type == "misc_feature":
            for qual in feature.qualifiers["note"]:
                if "CAGE" in qual:
                    found_features.append(feature)


for record in SeqIO.parse(template_gbk, "genbank"):
    for feature in found_features:
        record.features.append(feature)
    


    SeqIO.write(record, output_handle, "genbank")

output_handle.close()

