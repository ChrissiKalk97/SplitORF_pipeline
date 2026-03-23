import os
import sys
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# protein_seqs_fasta = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_26.01.2026-13.07.17_NMD_for_paper/Proteins_for_masspec.fa'
# transcript_type = 'NMD'
# outfile = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_26.01.2026-13.07.17_NMD_for_paper/Degrons/DegronMD_analysis_NMD.csv'

protein_seqs_fasta = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Input2023/TSL_filtered/protein_coding_peptide_sequences_110_tsl_refseq_filtered.fa'
transcript_type = 'protein_coding'
outfile = '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_26.01.2026-13.07.17_NMD_for_paper/Degrons/DegronMD_analysis_protein_coding.csv'

outdir = os.path.dirname(outfile)
os.makedirs(outdir, exist_ok=True)


degron_seq_df = pd.read_excel(
    '/Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/check_results/Degrons/Degron_supp_table_S2.xlsx', header=1)

degron_seq_df = degron_seq_df[degron_seq_df['Class'] == 'DEGRON']
protein_seqs = list(SeqIO.parse(protein_seqs_fasta, 'fasta'))

degron_summary = np.zeros((len(protein_seqs), 2))
degron_output_frame = pd.DataFrame(columns=[
                                   'tID', 'degron sequence', 'nr of occurrences', 'transcript length', 'perc degron seq'])

frame_counter = 0
counter = 0
rownames = []
for record in protein_seqs:
    for degron in degron_seq_df['Degron Hit'].unique():
        # degron_pat = re.compile(degron)
        # occurrences = degron_pat.findall(str(record.seq))
        occurrences = str(record.seq).count(degron)
        degron_summary[counter, 0] = occurrences + degron_summary[counter, 0]
        degron_summary[counter, 1] = occurrences * \
            len(degron)/len(record.seq) + degron_summary[counter, 1]
        if occurrences > 0:
            degron_output_frame.loc[frame_counter] = [record.id, degron,
                                                      occurrences, len(record.seq), occurrences*len(degron)/len(record.seq)]
            frame_counter += 1
    counter += 1
    rownames.append(record.id)

# plotting
sns.histplot(degron_summary[:, 1]).set(xlabel='degron percentage',
                                       title='Degron percentage of each SplitORF')
plt.savefig(os.path.join(
    outdir, f'Degron_percentage_{transcript_type}_matrix.png'))
plt.close()

fig, ax = plt.subplots()
sns.histplot(degron_summary[:, 1], ax=ax)\
    .set(xlabel='degron percentage',
         title='Degron percentage of each SplitORF')
ax.set_xlim(0, 0.2)
plt.savefig(os.path.join(
    outdir, f'Degron_percentage_{transcript_type}_matrix_zoom.png'))
plt.close()


sns.histplot(degron_summary[:, 0]) .set(xlabel='Nr of degrons',
                                        title='Nr of degrons per SplitORF')
plt.savefig(os.path.join(
    outdir, f'Nr_of_degrons_{transcript_type}_matrix.png'))
plt.close()


fig, ax = plt.subplots()
sns.histplot(degron_summary[:, 0], ax=ax)\
    .set(xlabel='Nr of degrons',
         title='Nr of degrons per SplitORF')
ax.set_xlim(0, 50)
plt.savefig(os.path.join(
    outdir, f'Nr_of_degrons_{transcript_type}_matrix_zoom.png'))
plt.close()

fig, ax = plt.subplots()
sns.histplot(degron_summary[:, 0], ax=ax)\
    .set(xlabel='Nr of degrons',
         title='Nr of degrons per SplitORF')
ax.set_xlim(0, 25)
plt.savefig(os.path.join(
    outdir, f'Nr_of_degrons_{transcript_type}_matrix_zoom_25.png'))
plt.close()


degron_output_frame.to_csv(outfile)
