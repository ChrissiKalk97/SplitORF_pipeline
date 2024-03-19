#This script calls the aligning scripts to align the Ribo-Seq reads to the transcriptome (all transcripts downloaded from Ensembl)
# and calculates the overlap with the previously determined unique regions (main pipeline) of RI and NMD
#The following are the main directories for the pipeline outputs and for the alignment outputs

#define the directories for ourput
outputBowtie="./Output/BOWTIE_transcriptome"
nmd="./Output/BOWTIE_transcriptome/NMD_transcriptome"
ri="./Output/BOWTIE_transcriptome/RI_transcriptome"
#define directories for input: results SplitORf pipeline
nmdmain="./Output/run_18.12.2023-10.14.23_new_NMD_longest_isoforms_with_UCSC"
rimain="./Output/run_18.12.2023-12.09.08_new_RI_longest_isoforms_with_UCSC"

#The following are the paths to the riboseq reads for each sample
#controls
OHMX20220060_001="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_001.fastq"
OHMX20220060_002="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_002.fastq"
OHMX20220060_003="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_003.fastq"

#treatments
OHMX20220060_004="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_004.fastq"
OHMX20220060_005="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_005.fastq"
OHMX20220060_006="/Users/christina/Documents/Riboseq_Siragusa/bam/OHMX20220060_006.fastq"

#arrays with sample names
sample_array=("$OHMX20220060_002" "$OHMX20220060_003" "$OHMX20220060_004" "$OHMX20220060_005" "$OHMX20220060_006")
sample_array_full=("$OHMX20220060_001" "$OHMX20220060_002" "$OHMX20220060_003" "$OHMX20220060_004" "$OHMX20220060_005" "$OHMX20220060_006")

#Create a Logfile for the alignments in the output directory
exec > >(tee -i $outputBowtie/AlignmentLogfile.txt)
exec 2>&1

#The following block calls the Bowtie_Align_transcriptomic script, which creates a BOWTIE index for the transcriptome and aligns
#Ribo-seq data against it before checking the overlap with the determined unique regions. 
#For further analysis a file with random regions from the 3' and 5' UTR with the same length distribution as the unique regions
#is also created and used to determine background overlap

###everything in between should not be quoted, just to be faster#############################################
echo "Starting alignment against transcripts"
source ./Bowtie_Align_transcriptomic.sh -i ./Input2023/cDNA_fasta_Ensembl_transcripts_111.fa 10 $outputBowtie/Transcriptomic_Bowtie_index $OHMX20220060_001 $nmd/001_NMD\
 $nmdmain/Unique_DNA_Regions_for_comparison.bed ./Input2023/ExonPositions_correct.bed\
 ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered.fa #-i ./Input2023/cDNA_fasta_Ensembl_transcripts.fa #this is verion 111

echo "=====...............	8%"

counter=2
for file in "${sample_array[@]}"
do
source ./Bowtie_Align_transcriptomic.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $file $nmd/00${counter}_NMD\
 $nmdmain/Unique_DNA_Regions_for_comparison.bed ./Input2023/ExonPositions_correct.bed\
 ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered.fa
counter=$((counter+1))
done 

echo "==========..........	50%"

counter=1
for file in "${sample_array_full[@]}"
do
source ./Bowtie_Align_transcriptomic.sh 10 $outputBowtie/Transcriptomic_Bowtie_index $file\
 $ri/00${counter}_RI $rimain/Unique_DNA_Regions_for_comparison.bed ./Input2023/ExonPositions_correct.bed\
  ./Input2023/Genomic_positions.bed ./Input2023/UTRs_merged_filtered.fa
counter=$((counter+1))
done
echo "====================	100%"

#This block corrects the variables for the rmd input and calls two RMD scripts, who give insight into the alignment of the Ribo-seq data
nmd=$nmd/
ri=$ri/
outputBowtie="./Output/BOWTIE_transcriptome"
R -e "rmarkdown::render('RiboSeqReportTranscriptomic.Rmd',output_file='./Output/BOWTIE_transcriptome/RiboSeq_Report.html',params=list(args = c('$outputBowtie', '$nmd','$ri')))"
#R -e "rmarkdown::render('RiboSeqReportGenomic.Rmd',output_file='./Output/RiboSeq_Report_genomic.html',params=list(args = c('$outputBowtie','$outputBowtie')))"