import re


vector = 'standard plasmid (grna only) Addgene# 6577 (CAGE SP78)'

test = re.search(r'\(C.*\)|\(P.*\)',vector)
print(test.group()) # (CAGE SP78)

#maketrans()
vector_clean = test.group().replace("CAGE","").translate(str.maketrans("()","  ")).strip()



print(vector_clean) # CAGE SP78



vector_file_dict = {
"SP78":	"SP78-AG65777-BPK1520 w LacZ insert_EMPTY.dna",
"SNP28": "SNP28-LentiPuro-AG52963 Empty.dna",
"PC292": "PC292_AG96925_pXPR_050_Lenti_empty.dna",
"PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
"PC291": "PC291_LentiPuro-P65-HSF_Empty.dna",
"SNP36": "SNP36 LentiCRISPRv2 -AG52961 Empty.dna",
"SS363": "SS363 pAW12.lenti.GFP _#104374_empty.dna",
"SNP112": "SNP112_AG82416_LentiCRISPRv2GFP Empty.dna",
"SP16":	"SP16 - LRG (Lenti_sgRNA_EFS_GFP) - AG65656_empty.dna",
"PC529": "PC529_LVA_plent-gRNA-Ametrine_Hongbo_lab_EMPTY_20Ns_update.dna",
"PC329": "PC329_pLentiCRIPSRv2_mCherry_AG99154_empty.dna",   
"PC226": "PC226-LMA_Hongbo_Empty.dna",
"SP76":	"SP76-LV04 U6gRNA PPB Sanger Crispr Vector_MP_newEMPTY_MP.dna",
"JK129": "JK129_AG61427_lenti sgRNA-MS2-Zeo_Empty.dna"
}


vector_file = vector_file_dict.get(vector_clean)

print(vector_file) # PC78-AG65777-BPK1520 w LacZ insert_EMPTY.dna