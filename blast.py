from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

result = NCBIWWW.qblast("blastn", "nt", "ATGGTATTAAAAGAGAATTTGCGCGTGAAGCAGGTACACCAGGGGGAGGACCCCATGGATATTGGTTACCAGCCGGGGAATGTTTGGTCCGTCGGCTCGGAGGTGGATGTGAGCTGTTGTGTATTTGTGGAACGAGCGGTGAAGCCTGCCGAAATTAGCGGAACCGTCAGAGAATGCCGCCACGCGCCTGCGGGAAACAGTTGCAACGCCGATTGTTTGTCGTCACCGACTGTGCACGTCATTATTTTCCTCACAACGACCACCGAGGGAGGCGTGTGCCAGATGGCGAAGGCGAAGTACGTTGGCCGCCGCCAAGAAGCGATCTGTCTGCAGGCCGCCCTCGGAGAAATGCTGAGCGAGATCACAAATCCGCTATGA")

# record = list(NCBIXML.parse(result))
# print(record)

with open("text.xml", "w+") as file:
    file.write(result.read())
    