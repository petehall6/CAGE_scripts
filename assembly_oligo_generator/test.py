import pandas as pd
from Bio.Seq import Seq


fwd_backbone = 'CACCG'
rev_backbone = 'AAAC'

my_seq = Seq('ATTATAAATACCGGCCCCGG')

rev_my_seq = my_seq.reverse_complement()

print(str(fwd_backbone + my_seq))
print(str(rev_backbone + rev_my_seq + 'C'))





