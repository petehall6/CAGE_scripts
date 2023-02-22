from Bio import Entrez,SeqIO
import os

os.chdir(os.path.dirname(os.path.abspath(__file__)))

Entrez.email = "pete.hall@stjude.org"

#get geneID
#get gene_summary
#fetch_gbk_record from summary

gene = input("Enter gene target gene: ").upper()

esearch_term = "{org}[Orgn] AND {gene}[Gene] AND alive[property]".format( 
    org = "HUMAN",
    gene=gene
)


def get_geneID():
    gbk_ent_handle = Entrez.esearch(db="gene", term=esearch_term,rettype="gb", retmode='text')
    gbk_ent_record = Entrez.read(gbk_ent_handle, "genbank")
    gene_id = gbk_ent_record["IdList"][0]

    return gene_id

def get_geneSummary(gene_id):
    
    rec = Entrez.read(Entrez.esummary(id=gene_id, db="gene"))
    
    return rec["DocumentSummarySet"]["DocumentSummary"][0]

def get_gbk(gene_id, extend_length=0):
    gene_id = int(gene_id)
    
    gene_summary = get_geneSummary(gene_id)

    genomic_info = gene_summary["GenomicInfo"]    
    
    genomic_rec_id = genomic_info[0]["ChrAccVer"]
    
    seq_start = int(genomic_info[0]["ChrStart"]) +1
    seq_stop = int(genomic_info[0]["ChrStop"]) +1
    
    strand = "1" if seq_start < seq_stop else "2"
    
    seq_start,seq_stop = sorted([seq_start, seq_stop])
    
    seq_start -= int(extend_length)
    seq_stop += int(extend_length)
    
    entrez_keywords = {
        "id": genomic_rec_id,
        "rettype": "gbwithparts",
        "retmode": "text",
        "seq_start": seq_start,
        "seq_stop": seq_stop,
        "strand": strand,
    }
    print("Fetching target genbank file.\n")
    gbk_handle = Entrez.efetch(db="nuccore", **entrez_keywords)
    
    gbk_rec = SeqIO.read(gbk_handle, format="genbank")
    gbk_handle.close()
    
    SeqIO.write([gbk_rec], open(f"{gene}_custom.gbk", "w"), "genbank")
    print("Completed")


gene_id = get_geneID()

get_geneSummary(gene_id)

get_gbk(gene_id, extend_length=0)





