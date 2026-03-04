#!/bin/bash


splitorf_pipeline="/Users/christina/Documents/bioconda-splitorf-pipeline"



#################################################################################
# ------------------ RUN SO PIPELINE                         ------------------ #
#################################################################################

eval "$(conda shell.bash hook)"
conda activate SplitORF
cd $splitorf_pipeline
bash $splitorf_pipeline/run_splitorfs_pipeline.sh \
 $splitorf_pipeline/Input2023/protein_coding_peptide_sequences_tsl_eq_filtered_29_09_25.fa \
 $splitorf_pipeline/Input2023/NMD_transcripts_CDNA.fa \
 $splitorf_pipeline/Input2023/ENSEMBLhuman110PFAMformat.bed \
 $splitorf_pipeline/Input2023/protein_coding_transcript_and_gene_cDNA_tsl_eq_filtered_29_09_25.fa \
 $splitorf_pipeline/Input2023/ExonCoordsWIthChr110.bed \
 diamond \
 $splitorf_pipeline/Input2023/Ens_110_CDS_coordinates_genomic_protein_coding_tsl_refseq_filtered.bed \
 $splitorf_pipeline/Output