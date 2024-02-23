# This script fetches the genomic positions of all transcripts and aplies them to the positions of the found unique regions
# It takes as first Input a tsv with the genomic transcript positions in the following format:
# Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp) strand
# The second Input file is the unique DNA regions file produced by the main pipeline
# The third argument needed for this script is the name of the output file with chromosome/scaffold name
# The fourth argument needed for this script is the name of the output file with geneID|transcriptID
import sys
file = open(sys.argv[1],'r')
list={}
for line in file:
    elems = line.split("\t")
    ids=elems[0]+"|"+elems[1]
    temp=[elems[2],elems[3],elems[4],elems[5]]
    list[ids]=temp
uniquefile = open(sys.argv[2],'r')
uniquelist=[]
for line in uniquefile:
    elements = line.strip()
    elements = elements.split("\t")
    orfIDs = elements[0].split(":")
    #start = int(orfIDs[2])
    #end = int(orfIDs[3])
    #temp=[orfIDs[0],start,end,orfIDs[1],orfIDs[2],orfIDs[3]]
    temp=[orfIDs[0],elements[1],elements[2],orfIDs[1]]
    uniquelist.append(temp)
with open(sys.argv[3],'w') as f:
    with open(sys.argv[4],'w') as f2:
        for i in uniquelist:
            if(int(list[i[0]][3])==1):
                i[1] = int(i[1]) + int(list[i[0]][1])
                i[2] = int(i[2]) + int(list[i[0]][1])
            else:
                st = i[1]
                en = i[2]
                i[2] = int(list[i[0]][2]) - int(st)
                i[1] = int(list[i[0]][2]) - int(en)
            f.write(list[i[0]][0] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n") #i[0] + "\t" +  set at start to get transcript and gene id
            #f2.write(i[0] + ":" + i[3] + ":" + i[4] + ":" + i[5] +"\t" + str(i[1]) + "\t" + str(i[2]) + "\n") #i[0] + "\t" +  set at start to get transcript and gene id
            f2.write(i[0] + ":" + i[3] + "\t" + str(i[1]) + "\t" + str(i[2]) + "\n")  # i[0] + "\t" +  set at start to get transcript and gene id

