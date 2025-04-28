import os
from Bio import Entrez, SeqIO

SIZE = 5000

def clean_files():
    
    holding_dir = os.path.join(os.getcwd(), "previous_gbks")
    
    for file in os.listdir(os.getcwd()):
        if file.endswith(".gbk"):
            try:
                os.rename(file, os.path.join(holding_dir, file))
            except FileNotFoundError:
                pass

def gene_ID_to_summary(geneID):
    """
    MAIN FAILURE POINT OF ENTIRE PROGRAM. NCBI XML STYLE
    CHANGES WILL CAUSE THIS TO START THROWING KEYERRORS / INDEXERRORS
    """
    geneID = int(geneID)
    assert geneID >= 0, "Negative gene IDs not allowed : %s" % geneID
    try:
        rec = Entrez.read(Entrez.esummary(id=geneID, db="gene"))
    except RuntimeError:
        raise Exception("NCBI gene ID %s does not exist" % geneID)

    # strip away the style of the raw record and just return the dictionary
    # that has the actual gene info
    return rec["DocumentSummarySet"]["DocumentSummary"][0]

print("Cleaning files...\n\n")

#clean_files()

gene_name = input("Enter gene name: ")
species= input("Enter species human or mouse: ")

Entrez.email = 'peter.hall@stjude.org'

if species == 'h':
    species = 'human'
elif species == 'm':
    species = 'mouse'

esearch_term = "{sp}[Orgn] AND {gene}[Gene] AND alive[property]".format(
    sp=species, gene=gene_name
)

gene_id_list = Entrez.read(
    Entrez.esearch("gene", term=esearch_term, retmax=20)
)["IdList"]


if len(gene_id_list) == 1:
    gene_id = gene_id_list[0]
else:
    gene_id = input(f"Select gene_id from list: {gene_id_list}\n")



"""
gene_id is an integer
SIZE is the number of bases to include upstream AND downstream
of the gene region as defined by ncbi
"""
gene_id = int(gene_id)

gene_summary = gene_ID_to_summary(gene_id)

genomic_info = gene_summary["GenomicInfo"]
assert len(genomic_info) == 1
# ID for the primary reference chromosome the gene is located on
genomic_record_id = genomic_info[0]["ChrAccVer"]

# for some reason these coords are different by one than what you get fromy
# manually searching ncbi gene database
seq_start = int(genomic_info[0]["ChrStart"]) + 1
seq_stop = int(genomic_info[0]["ChrStop"]) + 1

strand = "1" if seq_start < seq_stop else "2"

# Biopython freaks out if seq_stop < seq_start
seq_start, seq_stop = sorted([seq_start, seq_stop])

seq_start -= int(SIZE)
seq_stop += int(SIZE)

entrez_keywords = {
    "id": genomic_record_id,
    "rettype": "gbwithparts",
    "retmode": "text",
    "seq_start": seq_start,
    "seq_stop": seq_stop,
    "strand": strand,
}

for word in entrez_keywords:
    print(f"{word}: {entrez_keywords[word]}")

try:
    gbk_entrez_handle = Entrez.efetch(db="nuccore", **entrez_keywords)

    gbk_record = SeqIO.read(gbk_entrez_handle, format="genbank")
    gbk_entrez_handle.close()

    gbk_record

    SeqIO.write([gbk_record], open(f"gene_{gene_name}.gbk", "w"), "genbank")
    
    print(f"Gene written to gene_{gene_name}.gbk")
    print("\nDone")
    print("*" * 50)
except Exception as e:
    print(e)
    

