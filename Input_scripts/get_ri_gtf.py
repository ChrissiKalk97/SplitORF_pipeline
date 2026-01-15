#### filters a FASTA file for transcripts present in a gtf file#####
from pygtftk.gtf_interface import GTF
from Bio import SeqIO
import argparse
from typing import List


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    parser.add_argument("--ri_transcript_fasta",
                        help="Path to RI transcript FASTA")

    parser.add_argument("--ensembl_gtf_path",
                        help="Path to full Ensembl GTF")

    parser.add_argument("--outname",
                        help="Name of output GTF")

    return parser.parse_args()


def main(ri_transcript_fasta, ensembl_gtf_path, outname):
    def get_transcript_string(transcript_ids: List[str]) -> str:
        '''get string of transcript or other ids for filtering 
        of a GTF class object from pygtftk'''
        transcript_string = ''
        if len(transcript_ids) > 1:
            for transcript_id in transcript_ids[:-1]:
                transcript_string += transcript_id+','
        transcript_string += transcript_ids[-1]
        return transcript_string

    ri_transcripts = []
    for record in SeqIO.parse(ri_transcript_fasta, "fasta"):
        ri_transcripts.append(record.description.split('|')[1])

    ensembl_gtf = GTF(ensembl_gtf_path, check_ensembl_format=False)
    ri_transcript_string = get_transcript_string(ri_transcripts)

    ri_gtf = ensembl_gtf.select_by_key('feature', 'transcript,exon').select_by_key(
        'transcript_id', ri_transcript_string)

    ri_gtf.write(outname)


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    ri_transcript_fasta = args.ri_transcript_fasta
    ensembl_gtf_path = args.ensembl_gtf_path
    outname = args.outname

    main(ri_transcript_fasta, ensembl_gtf_path, outname)
