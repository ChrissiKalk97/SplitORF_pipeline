#!/bin/bash

# ----- This script is the master script to run the extended split-ORF pipeline. The Inputs are described below, ----- #
# ----- the outputs are all stored in a new folder which is created at the given location. The most interesting  ----- #
# ----- outputs are                                                                                              ----- #

# ----- Help message: ----- #
usage="
Usage: ./run_pipeline.sh [-h] proteins.fa transcripts.fa annotation.bed longestproteincodingtranscripts_DNA.fa longestproteincodingtranscripts_Protein.fa GenomicPositions.bed exonPositions.bed

proteins.fa 			should be a multi fasta file containing the amino acid sequences of the proteins that are used as reference (whole transcriptome).
transcripts.fa 			should be a multi fasta file containing the DNA sequences of the reads/transcripts that shall be analyzed.
annotation.bed 			should be a bedfile containing the annotations for the used genome build.
				The standard annotation files for human and mouse (ENSEMBL 95) can be found in the annotations directory in the SplitORF directory.
proteincodingtranscripts   CDS region as DNA sequence of all protein coding transcripts, can be downloaded from Ensmebl by selecting Sequence and coding regions
GenomicPositions.bed		A bed file containing the chromosome specific positions of all transcripts
should be in the format: Gene stable ID	Transcript stable ID	Chromosome/scaffold name	Transcript start (bp)	Transcript end (bp) Strand
exonPositions.bed			A bed file containing the chromosome specific position of all exons
should be in the format: Gene stable ID	Transcript stable ID	Exon region start (bp)	Exon region end (bp)	Transcript start (bp)	Transcript end (bp)	Strand
Tipp: The longest Isoforms can be extracted using the getLongestIsoform.py script within the Uniqueness_scripts folder (fastas must be sorted beforehand e.g. with seqkit sort)

where:
-h	show this help"

# ----- available options for the programm -----#
while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done

# ----- Check that all arguments are provided and are in the correct file format.                         ----- #
# ----- Give error messages according to errors in the call and exit the programm if something goes wrong ----- #
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages

# ----- check for right number of arguments ----- #
if [ "$#" -ne 6 ]; then 
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
  
# ----- check if any argument is a directory ----- #
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ] || [ -d "$5" ] || [ -d "$6" ] ; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
  
# ----- check if every argument has the correct file extension ----- #
elif ! [[ $1 =~ \.fa$ ]] || ! [[ $2 =~ \.fa$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.fa$ ]] || ! [[ $5 =~ \.bed$ ]]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi

# ----- Assign the given arguments to variables ----- #
proteins=$1
transcripts=$2
annotation=$3
proteinCodingTranscripts=$4
exonPositions=$5
align_method=$6

# ----- create the directory "Output" if it does not already exist ----- #
if [ -d "./Output" ]
then
	echo "Directory ./Output already exists."
else
	echo "Directory ./Output does not exist, creating ..."
	mkdir ./Output
	echo "Directory ./Output created."
fi

# ----- create a run specific output folder ----- #
timestamp=$(date "+%d.%m.%Y-%H.%M.%S")
if [ -d "./Output/run_$timestamp" ]
then 
	echo "Directory ./Output/run_$timestamp already exists." >&2
	exit 1
else
	mkdir ./Output/run_$timestamp
fi

# ----- Create a log file that saves all text outputs of the pipeline ----- #
exec > >(tee -i ./Output/run_$timestamp/Logfile.txt)
exec 2>&1

echo "Log Location should be: [ ./Output/run_$timestamp ]"

# ----- initialize conda ----- #
source $(conda info --base)/etc/profile.d/conda.sh

# ----- activate conda SplitORF environment needed for the pipeline ----- #
conda activate SplitORF
output="./Output/run_$timestamp"
echo "*********$output**********"
echo "run the Pipeline with: " $output $proteins $transcripts $annotation $proteinCodingTranscripts $exonPositions

# ----- create Orf sequences using the OrfFinder script ----- #
echo "searching for possible ORFs"
python ./SplitOrfs-master/OrfFinder_py3.py $transcripts > $output/OrfProteins.fa

# ----- Determine Unique Protein regions by calling mummer maxmatch with a minimum length of 10, annotating the ----- #
# ----- matches (non-unique regions) in a bedfile and using bedtools subtract to get the non matching regions   ----- #
# ----- which are then annotated as the unique regions in another bedfile                                       ----- #
echo "Align ORF-transcripts(Protein) to protein coding transcripts mummer -maxmatch -l 6"
#proteins instead of proteins2: match against all protein coding transcripts and not only the longest ones!!!
mummer -maxmatch -l 6 $proteins ./Output/run_$timestamp/OrfProteins.fa > ./Output/run_$timestamp/Proteins_maxmatch_l8.mums

echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/Proteins_maxmatch_l8.mums ./Output/run_$timestamp/Protein_non_unique.bed ./Output/run_$timestamp/OrfProteins.fa ./Output/run_$timestamp/OrfProteins.bed ./Output/run_$timestamp/Unique_Protein_Regions.bed
python ./Uniqueness_scripts/Filter_small_regions.py $output/Unique_Protein_Regions.bed 6 $output/Unique_Protein_Regions_gt8.bed

echo "Select only the non-unique regions for the upcoming blast"
bedtools subtract -a $output/OrfProteins.bed -b $output/Unique_Protein_Regions_gt8.bed > $output/ProteinRegionsforBlast.bed
bedtools getfasta -fi $output/OrfProteins.fa -bed $output/ProteinRegionsforBlast.bed -fo $output/ProteinsforBlast.fa

if [ "$align_method" == "blast" ]
then
    # ----- create a blast database for all proteins using BLAST+ ----- #
    echo "Create BlastDB using the protein coding peptide fasta"
    makeblastdb -in $proteins -out $output/ProteinDatabase -dbtype prot

    # ----- use BlastP to align the translated ORFs to the proteins (currently using 20 threads), ----- #
    # ----- change -num_threads otherwise. -outfmt "6std" results in file format:                 ----- #
    # ----- qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  ----- #
    echo "Blast the non-unique regions against the DB"
    blastp -query $output/ProteinsforBlast.fa -db $output/ProteinDatabase -out $output/OrfsAlign.txt -evalue 10 -num_threads 8 -outfmt "6 std"

    # ----- sort the blastp output by the second column to group by proteins and delete unsorted file ----- #
    sort -k2 $output/OrfsAlign.txt > $output/OrfsAlign_sorted.txt
    rm $output/OrfsAlign.txt
elif [ "$align_method" == "diamond" ]
then
    #----- create a blast database for all proteins using BLAST+ ----- #
    echo "Create BlastDB using the protein coding peptide fasta"
    diamond makedb --in $proteins -d $output/ProteinDatabase

    # ----- use BlastP to align the translated ORFs to the proteins (currently using 20 threads), ----- #
    # ----- change -num_threads otherwise. -outfmt "6std" results in file format:                 ----- #
    # ----- qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore  ----- #
    echo "Blast the non-unique regions against the DB"
    diamond blastp -q $output/ProteinsforBlast.fa -d $output/ProteinDatabase -o $output/OrfsAlign.tsv --evalue 10 --very-sensitive

    # ----- sort the blastp output by the second column to group by proteins and delete unsorted file ----- #
    sort -k2 $output/OrfsAlign.tsv > $output/OrfsAlign_sorted.txt
    rm $output/OrfsAlign.tsv
fi

# ----- run the detection script to identify the valid ORFs according to the parameters set in the detection file ----- #
echo "Select the valid ORFs"
python ./SplitOrfs-master/DetectValidSplitOrfMatches_py3.py $output/OrfsAlign_sorted.txt > $output/ValidProteinORFPairs.txt

# ----- sort file per Orf-transcript ID on column 3. Here it is important to omit the head while sorting          ----- #
cat $output/ValidProteinORFPairs.txt | awk 'NR<2{print ;next}{print | "sort -k3"}'  > $output/ValidProteinORFPairs_sortCol3.txt

# ----- select the best ORF matches and convert to bedfile ----- #
python ./SplitOrfs-master/getLongestOrfMatches_py3.py $output/ValidProteinORFPairs_sortCol3.txt > $output/UniqueProteinORFPairs.txt
python ./SplitOrfs-master/makeBed_py3.py $output/UniqueProteinORFPairs.txt > $output/UniqueProteinMatches.bed

# ----- Add the domain overlap to the selected ORFs        ----- #
echo "Add functional overlap (PFAM)"
#intersect results are empty for the new run...
bedtools intersect -a $output/UniqueProteinMatches.bed -b $annotation -wa -F 1  -wb > $output/intersectResults.txt
python ./SplitOrfs-master/addFunctionalOverlap_py3.py $output/UniqueProteinORFPairs.txt $output/intersectResults.txt > $output/UniqueProteinORFPairs_annotated.txt

# ----- create a report displaying the results of the original pipeline steps ----- #
echo "Create Split-ORF report"

##need to remove . in front of the reltive path to the transcripts file, otherwise this file cannot be found when the path is put together
transcripts_R="${transcripts:1}"
echo "$transcripts_R"
R -e "rmarkdown::render('Split-ORF_Report.Rmd',output_file='./Output/run_$timestamp/Split-ORF_Report.html',params=list(args = c('/Output/run_$timestamp/ValidProteinORFPairs_sortCol3.txt','/Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt','$transcripts_R')))"


# ----- extract the Valid ORF sequences from the given transcripts ----- #
echo 'Select Split-ORF DNA Sequences'
#now look back at all ORFS, bc until now the unique regions were omitted
#when valid (unique best match) -> take the whole ORF in the new file: Valid_ORF_Proteins.bed
python ./Uniqueness_scripts/SelectValidOrfSequences.py ./Output/run_$timestamp/UniqueProteinORFPairs.txt ./Output/run_$timestamp/OrfProteins.bed ./Output/run_$timestamp/Valid_ORF_Proteins.bed
#via transcript IDs: select the valid transcript DNA sequences by knowing valid IDs 
python ./Uniqueness_scripts/Select_validORF_DNA_sequences.py $transcripts ./Output/run_$timestamp/Valid_ORF_Proteins.bed ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa 

# ----- Determine Unique DNA regions by calling mummer maxmatch with a minimum length of 20, annotating the matches ----- #
# ----- (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions which are     ----- #
# ----- then annotated as the unique regions in another bedfile                                                     ----- #
echo "Align ORF-transcripts(DNA) to protein coding transcripts with mummer -maxmatch"
mummer -maxmatch $proteinCodingTranscripts ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa > ./Output/run_$timestamp/DNA_maxmatch.mums
echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py ./Output/run_$timestamp/DNA_maxmatch.mums ./Output/run_$timestamp/DNA_non_unique.bed ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa ./Output/run_$timestamp/ValidORF_DNA_Sequences.bed ./Output/run_$timestamp/Unique_DNA_Regions.bed

# ----- Merge the bedfile entries if the start and end positions for the same transcript only differ by the length  ----- #
# ----- parameter of MUMmer or less. For more details see Merge_Bedfile.py                                          ----- #
python ./Uniqueness_scripts/Filter_small_regions.py ./Output/run_$timestamp/Unique_DNA_Regions.bed 20 ./Output/run_$timestamp/Unique_DNA_Regions_gt20.bed

# ----- Extract the valid ORF-Proteins annotated in UniqueProteinORFPairs.txt                               ----- #
echo "Extract only the valid regions"
python ./Uniqueness_scripts/SelectValidOrfSequences.py ./Output/run_$timestamp/UniqueProteinORFPairs.txt ./Output/run_$timestamp/Unique_Protein_Regions_gt8.bed ./Output/run_$timestamp/Unique_Protein_Regions_gt8_valid.bed

# ----- Filter the protein and DNA unique regions for weird SplitORF positions
python ./Uniqueness_scripts/filter_unique_regions.py ./Output/run_$timestamp/UniqueProteinMatches.bed ./Output/run_$timestamp/Unique_Protein_Regions_gt8_valid.bed protein /Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_$timestamp
python ./Uniqueness_scripts/filter_unique_regions.py ./Output/run_$timestamp/UniqueProteinMatches.bed ./Output/run_$timestamp/Unique_DNA_Regions_gt20.bed DNA /Users/christina/Documents/ExtendedSplitORFPipeline-master_13_11_23/Output/run_$timestamp

# ----- Use bedtools getfasta to extract the fasta sequences of the unique regions annotated in the produced bedfiles ----- #
echo "Finishing Steps"
bedtools getfasta -fi ./Output/run_$timestamp/ValidORF_DNA_Sequences.fa -fo ./Output/run_$timestamp/Unique_DNA_regions.fa -bed ./Output/run_$timestamp/Unique_DNA_Regions_gt20_filtered.bed
bedtools getfasta -fi ./Output/run_$timestamp/ORFProteins.fa -fo ./Output/run_$timestamp/Unique_Protein_regions.fa -bed ./Output/run_$timestamp/Unique_Protein_Regions_gt8_valid_filtered.bed

# ----- Reorganize Unique_DNA_Regions_gt20.bed for later intersection with riboseq Alignment ----- #
#the position of the unique region is given with respect to the transcript and not the ORF!
python ./Uniqueness_scripts/Bedreorganize_adapted.py ./Output/run_$timestamp/Unique_DNA_Regions_gt20_filtered.bed ./Output/run_$timestamp/Unique_DNA_Regions_for_riboseq.bed
python ./Uniqueness_scripts/Bedreorganize_Proteins.py ./Output/run_$timestamp/Unique_Protein_Regions_gt8_valid_filtered.bed ./Output/run_$timestamp/Unique_Protein_Regions_transcript_coords.bed

exonPositionsSorted=$(basename $exonPositions .bed)_sorted.bed
awk '$7 == "1"' $exonPositions | sort -k1,1 -k2,2 -k3,3n > ./Output/run_$timestamp/plus_strand.bed
awk '$7 == "-1"' $exonPositions | sort -k1,1 -k2,2 -k4,4nr > ./Output/run_$timestamp/minus_strand.bed
cat ./Output/run_$timestamp/plus_strand.bed ./Output/run_$timestamp/minus_strand.bed > $exonPositionsSorted
#rm ./Output/run_$timestamp/plus_strand.bed
#rm ./Output/run_$timestamp/minus_strand.bed

# ----- Get the genomic positions for the unique DNA and protein regions ----- #
echo "Change the Positions of the unique DNA and Protein as well as the ValidORF bed to their positions within the unspliced transcript"
exonPositionsTranscriptPositions=$(basename $exonPositions .bed)_transcript_positions.bed
echo "Change the positions to their genomic equivalent and transform into UCSC format"
python ./Genomic_scripts_18_10_24/ExonToTranscriptPositions.py $exonPositionsSorted $exonPositionsTranscriptPositions
#rm $exonPositionsSorted
python ./Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py  ./Output/run_$timestamp/Unique_DNA_Regions_for_riboseq.bed $exonPositionsTranscriptPositions ./Output/run_$timestamp/Unique_DNA_regions_genomic.bed
python ./Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py  ./Output/run_$timestamp/Unique_Protein_Regions_transcript_coords.bed $exonPositionsTranscriptPositions ./Output/run_$timestamp/Unique_Protein_regions_genomic.bed


# ----- calculate overlap between unique DNA and protein regions genomic         ----- #
echo "calculate genomic overlap between unique DNA and protein regions"
bedtools intersect -a ./Output/run_$timestamp/Unique_DNA_regions_genomic.bed -b ./Output/run_$timestamp/Unique_Protein_regions_genomic.bed > ./Output/run_$timestamp/Unique_Regions_Overlap.bed
# ----- calculate overlap between unique DNA and protein regions transcriptomic         ----- #
echo "calculate transcriptomic overlap between unique DNA and protein regions"
bedtools intersect -a ./Output/run_$timestamp/Unique_DNA_Regions_for_riboseq.bed -b ./Output/run_$timestamp/Unique_Protein_Regions_transcript_coords.bed > ./Output/run_$timestamp/Unique_Regions_Overlap_transcriptomic.bed

# ----- Create a report with basics statistics of the uniqueness scripts                       ----- #
R -e "rmarkdown::render('Extended_Pipeline_new.Rmd',output_file='./Output/run_$timestamp/Uniqueness_Report.html',
params=list(args = c('/Output/run_$timestamp/Unique_DNA_Regions.fa', 
'/Output/run_$timestamp/Unique_Protein_Regions.fa',
'/Output/run_$timestamp/Unique_DNA_Regions_gt20_filtered.bed', 
'/Output/run_$timestamp/Unique_Protein_Regions_gt8_valid_filtered.bed',
'/Output/run_$timestamp/UniqueProteinORFPairs.txt')))"



