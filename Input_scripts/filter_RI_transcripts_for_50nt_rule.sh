#!/bin/bash


# =============================================================================
# Script Name: 
# Description: filters RI transcripts as Input to the SO pipeline for 50nt rule
# Usage:       bash 
# Author:      Christina Kalk
# Date:        2026-01-13
# =============================================================================


# ----- initialize conda ----- #
source $(conda info --base)/etc/profile.d/conda.sh
conda activate pygtftk


so_input_dir="/home/ckalk/tools/SplitORF_pipeline/Input2023"

ri_transcripts="${so_input_dir}/RI_transcripts_CDNA.fa"
ensembl_full_gtf="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.113.chr.gtf"
ensembl_filtered_gtf="/projects/splitorfs/work/reference_files/filtered_Ens_reference_correct_29_09_25/Ensembl_110_filtered_equality_and_tsl1_2_correct_29_09_25.gtf"
genome_fasta="/projects/splitorfs/work/reference_files/Homo_sapiens.GRCh38.dna.primary_assembly_110.fa"
fiftynt_pipeline="/home/ckalk/tools/NMD_fetaure_composition"

ri_gtf="${so_input_dir}/RI_110.gtf"

cd /home/ckalk/tools/SplitORF_pipeline/Input_scripts
if [ ! -e "${ri_gtf}" ]; then
    python get_ri_gtf.py \
    --ri_transcript_fasta "${ri_transcripts}" \
    --ensembl_gtf_path "${ensembl_full_gtf}" \
    --outname "${ri_gtf}"
fi


cd "${fiftynt_pipeline}"
if [ ! -e "${fiftynt_pipeline}/Output/RI_transcripts_Ens" ]; then
    python fifty_nt_rule.py \
        "${ri_gtf}" \
        "${ensembl_filtered_gtf}" \
        "${genome_fasta}" \
        RI_transcripts_Ens_110
fi


cd -
python filter_RI_transcripts_for_50nt.py \
 --ri_transcript_fasta "${ri_transcripts}" \
 --fiftynt_csv "${fiftynt_pipeline}/Output/RI_transcripts_Ens/RI_transcripts_Ens.csv" \
 --outname "${so_input_dir}/RI_transcripts_50nt_CDNA.fa"