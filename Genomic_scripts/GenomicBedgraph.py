# This script fetches the genomic positions of all transcripts and aplies them to the positions of the found unique regions
# It takes as first Input a tsv with the genomic transcript positions in the following format:
# Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp)
# The second Input file is the bedgraph file
# The third argument needed for this script is the name of the output file
import sys
file = open(sys.argv[1],'r')
list={}
for line in file:
    elems = line.split("\t")
    ids=elems[0]+"|"+elems[1]
    temp=[elems[2],elems[3],elems[4],elems[5]]
    list[ids]=temp
graphfile = open(sys.argv[2],'r')
graphlist=[]
for line in graphfile:
    elements = line.strip()
    elements = elements.split("\t")
    id = elements[0]
    start = elements[1]
    end = elements[2]
    cov = elements[3]
    temp=[id,start,end,cov]
    graphlist.append(temp)
with open(sys.argv[3],'w') as f:
    for i in graphlist:
        if (int(list[i[0]][3]) == 1):
            i[1] = int(i[1]) + int(list[i[0]][1])
            i[2] = int(i[2]) + int(list[i[0]][1])
        else:
            st = i[1]
            en = i[2]
            i[2] = int(list[i[0]][2]) - int(st)
            i[1] = int(list[i[0]][2]) - int(en)
        f.write(list[i[0]][0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\t" + str(i[3]) + "\n") #i[0] + "\t" +  set at start to get transcript and gene id

