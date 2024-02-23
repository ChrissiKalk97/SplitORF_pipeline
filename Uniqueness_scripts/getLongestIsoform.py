from Bio import SeqIO
import sys


source = SeqIO.parse(sys.argv[1], "fasta")
with open(sys.argv[2], "w") as out:
    previousGeneID = "Dummy"
    longestSeqLen = 0
    sequences = []
    i = 0
    for line in source:
        geneID = line.id.split("|")
        if (geneID[0] == previousGeneID):
            if (len(line.seq) > longestSeqLen):
                longestSeqLen = len(line.seq)
                sequences[i-1]=line
        else:
            longestSeqLen = len(line.seq)
            sequences.append(line)
            previousGeneID = geneID[0]
            i = i+1
    SeqIO.write(sequences, out, "fasta")