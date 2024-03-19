##############################################################################
###BackgroundRegions_adapted.py takes a bedfile and a fasta as input and
###generates random regions of the fasta file with the length distribution as 
###given in the bed file and returns them in bed format
##############################################################################

#usage: python BackgroundRegions_adapted.py in.bed in.fasta out.bed
from Bio import SeqIO
import sys
import random


def get_length_dist(bed_file):
    #read in bedfile and obtain length distribution
    lengthdistribution=[]
    for line in bed_file:
        elems = line.split('\t')
        temp=int(elems[2])-int(elems[1])
        lengthdistribution.append(temp)
    return lengthdistribution

def get_seq_dictionary(fasta):
    #read in fasta file and get keys, add a number in case of
        #duplicate keys
    random_seqs={}
    i = 0
    for line in fasta:
        random_seqs[line.id + '_' + str(i)] = line.seq
        i += 1
    return random_seqs

def get_random_fasta_ids(lengthdistribution, random_seqs):
    #select random keys for the fasta file, as many as there are 
        #entries in the bed file
    #do not allow for duplicate keys
    randomlist=[]
    for i in range(len(lengthdistribution)):
        random_key = random.choices(list(random_seqs.keys()), k=1)
        while random_key in randomlist:
            random_key = random.choices(list(random_seqs.keys()), k=1)
        randomlist.append(random_key[0])
    return randomlist


def write_bed_output(outname, randomlist, random_seqs, lengthdistribution):
    #for each length from the length distribution
        #get a random sample of the same length from the fasta sequence
    with open(outname, 'w') as out:
        i = len(randomlist) - 1
        while i > 0:
            length=lengthdistribution.pop()
            random_key = randomlist.pop()
            random_seq = random_seqs[random_key]
            if length < len(random_seq):
                start = random.choice(range(len(random_seq)-length))
                end = start + length
                out.write(random_key.split('_')[0] + '\t' + str(start) + '\t' + str(end) + '\n')
                i-=1 
            else:
                #if the sampled fasta sequence is too short, sample a new random fasta entry
                #and reappend the length to resample
                random_key = random.choices(list(random_seqs.keys()), k=1)
                while random_key in randomlist:
                    random_key = random.choices(list(random_seqs.keys()), k=1)
                randomlist.append(random_key[0])
                lengthdistribution.append(length)



if len(sys.argv) < 3:
        print('usage python BackgroundRegions_adapted.py in.bed in.fasta out.bed')
else :
    file = open(sys.argv[1],'r')
    lengthdistribution = get_length_dist(file)  
    
    fasta = SeqIO.parse(sys.argv[2], 'fasta')
    random_seqs =  get_seq_dictionary(fasta)

    randomlist = get_random_fasta_ids(lengthdistribution, random_seqs)
    
    write_bed_output(sys.argv[3], randomlist, random_seqs, lengthdistribution)
    
        