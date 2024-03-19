#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs
nmd="./Output/Bowtie_transcriptome_pooled_no_5_adapted_110/NMD_transcriptome"
nmdmain="./Output/check_results/positions_of_unique_Regions/weird_unique_positions/filtered_unique_regions"
ri="./Output/Bowtie_transcriptome_pooled_no_5_adapted_110/RI_transcriptome"
rimain="./Output/check_results/positions_of_unique_Regions/weird_unique_positions/filtered_unique_regions"

outputBowtie="./Output/Bowtie_transcriptome_pooled_no_5_adapted_110"
#The following are the paths to the riboseq reads for each sample
#controls
control="/Users/christina/Documents/Riboseq_Siragusa/bam/control.fastq"

#treatments
treatment="/Users/christina/Documents/Riboseq_Siragusa/bam/treatment_no_5.fastq"


#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputBowtie/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align script, which creates a BOWTIE index for thetranscriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against transcripts"
source ./Bowtie_Align_transcriptomic_adapted.sh -i ./Input2023/cDNA_fasta_Ensembl_transcripts_110.fa 10 $outputBowtie/Transcriptomic_Bowtie_index $control $nmd/control_NMD\
 $nmdmain/Unique_DNA_Regions_for_comparison_filtered_NMD.bed ./Input2023/ExonPositions_correct.bed\
 ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered_110.fa

echo "=====...............	25%"

source ./Bowtie_Align_transcriptomic_adapted.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $treatment $nmd/treatment_NMD\
 $nmdmain/Unique_DNA_Regions_for_comparison_filtered_NMD.bed ./Input2023/ExonPositions_correct.bed\
 ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered_110.fa



echo "==========..........	50%"


source ./Bowtie_Align_transcriptomic_adapted.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
 $ri/control_RI $rimain/Unique_DNA_Regions_for_comparison_filtered_RI.bed ./Input2023/ExonPositions_correct.bed\
  ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered_110.fa

source ./Bowtie_Align_transcriptomic_adapted.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $control\
 $ri/treatment_RI $rimain/Unique_DNA_Regions_for_comparison_filtered_RI.bed ./Input2023/ExonPositions_correct.bed\
  ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered_110.fa

echo "===================	100%"

#This block corrects the variables for the rmd input and calls two RMD scripts, who give insight into the alignment of the Ribo-seq data
nmd=$nmd/
ri=$ri/
outputBowtie="./Output/Bowtie_transcriptome_pooled_no_5_adapted_110/"
R -e "rmarkdown::render('RiboSeqReportTranscriptomic.Rmd',output_file='./Output/Bowtie_transcriptome_pooled_no_5_adapted_110/RiboSeqReportTranscriptomic.pdf',params=list(args = c('$outputBowtie', '$nmd','$ri')) ,encoding="UTF-8")"
#R -e "rmarkdown::render('./Output/BOWTIE_transcriptome_filtered/RiboSeqReport.Rmd',output_file='./Output/RiboSeq_Report.html',params=list(args = c('$nmd','$ri')))"
#R -e "rmarkdown::render('./Output/BOWTIE_transcriptome_filtered/RiboSeqReportGenomic.Rmd',output_file='./Output/RiboSeq_Report_genomic.html',params=list(args = c('$outputBowtie','$outputBowtie')))"