from Bio import SeqIO, SearchIO
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd

# Load Toxoplasma protein sequences from a FASTA file
toxo_proteins = list(SeqIO.parse('toxo.fasta', 'fasta'))

# Load Human protein sequences from a FASTA file
human_proteins = list(SeqIO.parse('human.fasta', 'fasta'))

# Function to run BLAST and parse results
def blast_and_parse(sequence, database="nr", program="blastp", e_value_thresh=0.01):
    # Run BLAST
    result_handle = NCBIWWW.qblast(program, database, sequence.seq)
    
    # Parse BLAST output
    blast_records = NCBIXML.read(result_handle)

    print(blast_records)
    
    # Extract results
    results = []
    for alignment in blast_records.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < e_value_thresh:
                result = {
                    "Toxo_Protein": sequence.id,
                    "Human_Protein": alignment.title.split()[1], # Extracting the protein ID from title
                    "E-value": hsp.expect,
                    "Identity_Percentage": (hsp.identities / hsp.align_length) * 100
                }
                results.append(result)
    
    return results

# Example: Running BLAST for the first Toxoplasma protein
example_results = blast_and_parse(toxo_proteins[0])

# Display the example results
print(example_results)
