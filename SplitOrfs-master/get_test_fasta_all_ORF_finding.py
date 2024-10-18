from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord



sequences = []
seq1 = Seq('ATGAAAATGTTTATGAAATAAAAATGAAAATGAAAATGAATGAAATGAAATGAATGA')
seq1 = SeqRecord(seq=seq1, id='TEST1', name='<unknown name>', description='<unknown description>', dbxrefs=[])
seq2 = Seq('ATGAAAATGTTTTGA')
seq2 = SeqRecord(seq=seq2, id='TEST2', name='<unknown name>', description='<unknown description>', dbxrefs=[])
sequences.append(seq1)
sequences.append(seq2)
with open("test_all_ORF_finder.fa", "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")