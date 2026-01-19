#### filters a FASTA file for transcripts present in a gtf file#####
from Bio import SeqIO
import pandas as pd
import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="."
    )

    parser.add_argument("--ri_transcript_fasta",
                        help="Path to RI transcript FASTA")

    parser.add_argument("--fiftynt_csv",
                        help="Path to the CSV results from the 50nt run")

    parser.add_argument("--outname",
                        help="Name of output FASTA file")

    return parser.parse_args()


def main(ri_transcript_fasta, fiftynt_csv, outname):

    fiftynt_csv_df = pd.read_csv(fiftynt_csv)
    fiftynt_csv_df['50_nt'] = fiftynt_csv_df['50_nt'].astype(int)
    fiftynt_csv_df_nmd = fiftynt_csv_df[fiftynt_csv_df['50_nt'] == 1]

    ri_transcripts_nmd = []
    for record in SeqIO.parse(ri_transcript_fasta, "fasta"):
        if record.description.split('|')[1] in fiftynt_csv_df_nmd['Unnamed: 0'].to_list():
            ri_transcripts_nmd.append(record)

    with open(outname, "w") as output_handle:
        SeqIO.write(ri_transcripts_nmd, output_handle, "fasta")


if __name__ == "__main__":
    # ------------------ CONSTANTS ------------------ #
    args = parse_args()

    ri_transcript_fasta = args.ri_transcript_fasta
    fiftynt_csv = args.fiftynt_csv
    outname = args.outname

    main(ri_transcript_fasta, fiftynt_csv, outname)
