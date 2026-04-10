#!/bin/bash -l
# this script runs the downstream analysis steps for the validated regions per sample
# teh NMd inhibited samples are treated separately, for these validation only if 2 ORFs are Ribo-covered
# 1. get all validated regions union and statistics for NMD and RI
# 2. check for RBP validation from RBPDB: the database needs to be downloaded and paths 
# indicated for this to run


eval "$(conda shell.bash hook)"
conda activate Riboseq

rbpdb=false
pfam=false

outdir="/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome/..."
if [[ ! -d "${outdir}/NMD" ]]; then
    python so_categorization_coverage_pipeline.py \
        --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
        --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
        --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome" \
        --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
        --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome" \
        --region_type "NMD" \
        --sample_type "control"

        python so_categorization_coverage_pipeline.py \
        --so_results "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/UniqueProteinORFPairs.txt" \
        --so_categorization_df "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/so_categorization_df.csv" \
        --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome" \
        --ur_path "/projects/splitorfs/work/split-orf-prediction/Output/run_07.04.2026-16.05.28_NMD_cont_subtraction/Unique_DNA_Regions_genomic_final.bed" \
        --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/test_Ribo_val_conda/NMD_genome" \
        --region_type "NMD" \
        --sample_type "NMD_inhibition"
fi

# if [[ ! -d "${outdir}/RI" ]]; then
#     python so_categorization_coverage_pipeline.py \
#         --so_results "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/UniqueProteinORFPairs.txt" \
#         --ribo_coverage_path "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/RI_genome" \
#         --ur_path "/home/ckalk/tools/SplitORF_pipeline/Output/run_26.01.2026-14.09.13_RI_for_paper/Unique_DNA_Regions_genomic_final.bed" \
#         --result_dir "/projects/splitorfs/work/Riboseq/Output/Riboseq_genomic_single_samples/resample_q10_expression_filter/SO_validated_set_analysis" \
#         --region_type "RI"
# fi


if [[ "${rbpdb}" == true ]]; then 
    # bash ../downstream_analysis_validated_URs/run_rbpdb_analysis.sh > ../downstream_analysis_validated_URs/outreports_of_runs/run_rbpdb_analysis.out 2>&1
    bash ../downstream_analysis_validated_URs/run_rbp2go_analysis.sh > ../downstream_analysis_validated_URs/outreports_of_runs/run_rbp2go_analysis.out 2>&1
fi

if [[ "${pfam}" == true ]]; then 
    bash ../downstream_analysis_validated_URs/check_validated_pfam_domains.sh > ../downstream_analysis_validated_URs/outreports_of_runs/check_validated_pfam_domains.out 2>&1
fi