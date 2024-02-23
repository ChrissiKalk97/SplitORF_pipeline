from Bio import SeqIO
import sys
import random

file = open(sys.argv[1],'r')
lengthdistribution=[]
for line in file:
    elems = line.split("\t")
    temp=int(elems[2])-int(elems[1])
    lengthdistribution.append(temp)

source = SeqIO.parse(sys.argv[2], "fasta")
list={}
for line in source:
    list[line.id]=line.seq

randomlist=[]
for i in range(len(lengthdistribution)):
    randomlist.append(random.choice(range(len(list.keys()))))
    #randomly select one of the chromosomes or the NMD/RI transcripts by numbering them

source = SeqIO.parse(sys.argv[2], "fasta")
with open(sys.argv[3], "w") as out:
    i=0
    for line in source:
        while i in randomlist:
            length=lengthdistribution.pop()
            if length < len(line.seq):
                start = random.choice(range(len(line.seq)-length))
                end = start + length
                out.write(line.id + "\t" + str(start) + "\t" + str(end) + "\n")
            else:
                out.write(line.id + "\t" + "0" + "\t" + str(len(line.seq)) + "\n")
            randomlist.remove(i)
        i+=1