import sys
import re
from Bio import SeqIO

# this script extracts the protein sequences of the valid and unique SplitORF transcripts
# for MassSpec search
# input files are:
# 1: file with unique protein regions (fasta)
# 2: file with all possible ORFs (fasta)
# 3: Output file for protein sequences (fasta)


def extract_Split_ORF_ids(file1):
    Unique_regions = open(file1, 'r')
    Split_orf_ids = {}
    for unique_region in SeqIO.parse(Unique_regions, 'fasta'):
        Split_ORF_info = re.split(r'[|:]', unique_region.id)
        if Split_ORF_info[1] in Split_orf_ids.keys():
            if Split_ORF_info[2] not in Split_orf_ids[Split_ORF_info[1]]:
                Split_orf_ids[Split_ORF_info[1]].append(Split_ORF_info[2])
        else:
            Split_orf_ids[Split_ORF_info[1]] = [Split_ORF_info[2]]
    Unique_regions.close()
    return Split_orf_ids


def get_fasta_sequences(Split_orf_ids, fasta_name, outname):
    fin = open(fasta_name, 'r')
    fout = open(outname, 'w')
    for record in SeqIO.parse(fin, 'fasta'):
        header_info = re.split(r'[|:]', record.id)
        if header_info[1] in Split_orf_ids.keys():
            for transcript_id, list_ORFs in Split_orf_ids.items():
                if list_ORFs:
                    if transcript_id == header_info[1]:
                        # print(transcript_id,  header_info[1])
                        for ORF in list_ORFs:
                            if ORF == header_info[2]:
                                # print(transcript_id,  header_info[1], ORF, header_info[2])
                                record.id = 'sp|' + record.id + '|'
                                SeqIO.write(record, fout, "fasta")
    fin.close()
    fout.close()


if len(sys.argv) < 4:
    print("usage comapre_gtf.py Unique_protein_regions.fa OrfProteins.fa outputfile.fa")
else:
    Split_ORF_ids = extract_Split_ORF_ids(sys.argv[1])
    get_fasta_sequences(Split_ORF_ids, sys.argv[2], sys.argv[3])
