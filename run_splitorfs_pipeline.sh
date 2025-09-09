#!/bin/bash

# ----- This script is the master script to run the extended split-ORF pipeline. The Inputs are described below, ----- #
# ----- the outputs are all stored in a new folder which is created at the given location. The most interesting  ----- #
# ----- outputs are                                                                                              ----- #

# ----- Help message: ----- #
usage="
Usage: ./run_pipeline.sh [-h] proteins.fa transcripts.fa annotation.bed \\
    longestproteincodingtranscripts_DNA.fa longestproteincodingtranscripts_Protein.fa \\
    GenomicPositions.bed exonPositions.bed

Arguments:
  proteins.fa
      A multi-FASTA file containing the amino acid sequences of the proteins
      that are used as references (whole transcriptome, downloaded from Ensembl version 110).

  transcripts.fa
      A multi-FASTA file containing the DNA sequences of the reads/transcripts
      that shall be analyzed (input transcripts).

  PFAMannotation.bed
      A BED file containing the annotations for the used genome build.
      These were downloaded from Ensembl and modified with a custom script.

  reference_transcripts.fa
      A multi-FASTA file of all reference transcripts (should not contain the input
      transcripts).

  exonPositions.bed
      A BED file containing the chromosome-specific positions of all exons.
      These can be downloaded from Biomart using the structures functionality.
      Format:
        Gene stable ID	Transcript stable ID	Exon region start (bp)	
        Exon region end (bp)	Transcript start (bp)	Transcript end (bp)
        	Strand	Chromosome/scaffold name
  alignmethod
      blast or diamnond: diamond is faster while blast recovers more hits


Options:
  -h    Show this help message.
"


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
elif ! [[ $1 =~ \.(fa|fasta)$ ]] || ! [[ $2 =~ \.(fa|fasta)$ ]] || ! [[ $3 =~ \.(txt|bed)$ ]] || ! [[ $4 =~ \.(fa|fasta)$ ]] || ! [[ $5 =~ \.(txt|bed)$ ]]; then
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
  mkdir ./Output/run_$timestamp/tests
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
echo "Align ORF-transcripts(Protein) to protein coding transcripts mummer -maxmatch -l 8"
#proteins instead of proteins2: match against all protein coding transcripts and not only the longest ones!!!
mummer -maxmatch -l 8 $proteins $output/OrfProteins.fa > $output/Proteins_maxmatch_l8.mums

echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py\
 $output/Proteins_maxmatch_l8.mums\
 $output/Protein_non_unique.bed\
 $output/OrfProteins.fa\
 $output/OrfProteins.bed\
 $output/Unique_Protein_Regions.bed

python ./Uniqueness_scripts/Filter_small_regions.py\
 $output/Unique_Protein_Regions.bed\
 8\
 $output/Unique_Protein_Regions_gt8.bed

echo "Select only the non-unique regions for the upcoming blast"
bedtools subtract -a $output/OrfProteins.bed -b $output/Unique_Protein_Regions_gt8.bed > $output/ProteinRegionsforBlast.bed
bedtools getfasta -fi $output/OrfProteins.fa -bed $output/ProteinRegionsforBlast.bed -fo $output/ProteinsforBlast.fa

if [[ "$align_method" == "blast" ]]; then
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
elif [[ "$align_method" == "diamond" ]]; then
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
python ./SplitOrfs-master/getLongestOrfMatches_py3.py\
 $output/ValidProteinORFPairs_sortCol3.txt\
 > $output/UniqueProteinORFPairs.txt

python ./SplitOrfs-master/makeBed_py3.py\
 $output/UniqueProteinORFPairs.txt\
 > $output/UniqueProteinMatches.bed

# ----- Add the domain overlap to the selected ORFs        ----- #
echo "Add functional overlap (PFAM)"
#intersect results are empty for the new run...
bedtools intersect -a $output/UniqueProteinMatches.bed -b $annotation -wa -F 1  -wb > $output/intersectResults.txt

python ./SplitOrfs-master/addFunctionalOverlap_py3.py\
 $output/UniqueProteinORFPairs.txt\
 $output/intersectResults.txt\
 > $output/UniqueProteinORFPairs_annotated.txt

# ----- create a report displaying the results of the original pipeline steps ----- #
echo "Create Split-ORF report"

##need to remove . in front of the reltive path to the transcripts file, otherwise this file cannot be found when the path is put together
if [[ $transcripts == ./* ]]; then
  transcripts_R="${transcripts:1}"
else
  transcripts_R="${transcripts}"
fi


R -e "rmarkdown::render('Split-ORF_Report.Rmd',
output_file='$output/Split-ORF_Report.html',
params=list(args = c('/Output/run_$timestamp/ValidProteinORFPairs_sortCol3.txt',
'/Output/run_$timestamp/UniqueProteinORFPairs_annotated.txt',
'$transcripts_R')))"


# ----- extract the Valid ORF sequences from the given transcripts ----- #
echo 'Select Split-ORF DNA Sequences'
# now look back at all ORFS, bc until now the unique regions were omitted
# when valid (unique best match) -> take the whole ORF in the new file: Valid_ORF_Proteins.bed
python ./Uniqueness_scripts/SelectValidOrfSequences.py\
 $output/UniqueProteinORFPairs.txt\
 $output/OrfProteins.bed\
 $output/Valid_ORF_Proteins.bed

# via transcript IDs: select the valid transcript DNA sequences by knowing valid IDs 
python ./Uniqueness_scripts/Select_validORF_DNA_sequences.py\
 $transcripts $output/Valid_ORF_Proteins.bed\
 $output/ValidORF_DNA_Sequences.fa 

# ----- Determine Unique DNA regions by calling mummer maxmatch with a minimum length of 20, annotating the matches ----- #
# ----- (non unique regions) in a bedfile and using bedtools subtract to get the non matching regions which are     ----- #
# ----- then annotated as the unique regions in another bedfile                                                     ----- #
# ----- 20 is the defualt parameter for -l matchlength                                                              ----- #
echo "Align ORF-transcripts(DNA) to protein coding transcripts with mummer -maxmatch"
mummer -maxmatch $proteinCodingTranscripts $output/ValidORF_DNA_Sequences.fa > $output/DNA_maxmatch.mums

echo "Select the non matching regions as unique regions and save as bedfile"
python ./Uniqueness_scripts/Find_Unique_Regions.py\
 $output/DNA_maxmatch.mums\
 $output/DNA_non_unique.bed\
 $output/ValidORF_DNA_Sequences.fa\
 $output/ValidORF_DNA_Sequences.bed\
 $output/Unique_DNA_Regions.bed

# ----- Merge the bedfile entries if the start and end positions for the same transcript only differ by the length  ----- #
# ----- parameter of MUMmer or less. For more details see Merge_Bedfile.py                                          ----- #
python ./Uniqueness_scripts/Filter_small_regions.py\
 $output/Unique_DNA_Regions.bed\
 20\
 $output/Unique_DNA_Regions_gt20.bed

# ----- Extract the valid ORF-Proteins annotated in UniqueProteinORFPairs.txt                               ----- #
echo "Extract only the valid regions"
python ./Uniqueness_scripts/SelectValidOrfSequences.py\
 $output/UniqueProteinORFPairs.txt\
 $output/Unique_Protein_Regions_gt8.bed\
 $output/Unique_Protein_Regions_gt8_valid.bed

# ----- Filter the protein and DNA unique regions for weird SplitORF positions
python ./Uniqueness_scripts/filter_unique_regions.py\
 $output/UniqueProteinMatches.bed\
 $output/Unique_Protein_Regions_gt8_valid.bed\
 protein\
 $output

python ./Uniqueness_scripts/filter_unique_regions.py\
 $output/UniqueProteinMatches.bed\
 $output/Unique_DNA_Regions_gt20.bed\
 DNA\
 $output

# ----- Use bedtools getfasta to extract the fasta sequences of the unique regions annotated in the produced bedfiles ----- #
echo "Finishing Steps"
bedtools getfasta\
 -fi $output/ValidORF_DNA_Sequences.fa\
 -fo $output/Unique_DNA_Regions.fa\
 -bed $output/Unique_DNA_Regions_gt20_filtered.bed

bedtools getfasta\
 -fi $output/ORFProteins.fa\
 -fo $output/Unique_Protein_Regions.fa\
 -bed $output/Unique_Protein_Regions_gt8_valid_filtered.bed

# ----- get the ORF seqeunces of those predicted proteins with unique regions  ----- #
python ./Uniqueness_scripts/get_protein_seqs_for_Masspec.py\
 $output/UniqueProteinORFPairs.txt\
 $output/OrfProteins.fa\
 $output/Proteins_for_masspec.fa

# ----- Reorganize Unique_DNA_Regions_gt20.bed for later intersection with riboseq Alignment ----- #
# Note: The position of the unique region is given with respect to the transcript and not the ORF!
python ./Uniqueness_scripts/Bedreorganize_adapted.py \
    $output/Unique_DNA_Regions_gt20_filtered.bed \
    $output/Unique_DNA_Regions_for_riboseq.bed

python ./Uniqueness_scripts/Bedreorganize_Proteins.py \
    $output/Unique_Protein_Regions_gt8_valid_filtered.bed \
    $output/Unique_Protein_Regions_transcript_coords.bed


exonPositionsSorted=$(basename $exonPositions .bed)_sorted.bed
awk '$7 == "1" || $7 == "+"' $exonPositions | sort -k1,1 -k2,2 -k3,3n > $output/plus_strand.bed
awk '$7 == "-1" || $7 == "-"' $exonPositions | sort -k1,1 -k2,2 -k4,4nr > $output/minus_strand.bed
cat $output/plus_strand.bed $output/minus_strand.bed > $output/$exonPositionsSorted
rm $output/plus_strand.bed
rm $output/minus_strand.bed


# ----- Get the genomic positions for the unique DNA and protein regions ----- #
echo "Change the Positions of the unique DNA and Protein as well as the ValidORF bed to their positions within the unspliced transcript"
exonPositionsTranscriptPositions=$(basename $exonPositions .bed)_transcript_positions.bed
echo "Change the positions to their genomic equivalent and transform into UCSC format"
python ./Genomic_scripts_18_10_24/ExonToTranscriptPositions.py\
 $output/$exonPositionsSorted\
  $output/$exonPositionsTranscriptPositions




# ----- Test functionality of the exon coordinate conversion ----- #
mkdir $output/tests
python ./Genomic_scripts_18_10_24/ExonToTranscriptPositions.py\
 ./Genomic_scripts_18_10_24/test/exon_transcript_positions_unit_test.bed\
  $output/tests/exon_transcript_positions_results.bed

python ./Genomic_scripts_18_10_24/test/test_transcript_exon_positions.py\
 $output/tests/exon_transcript_positions_results.bed\
 $output/$exonPositionsSorted\
 $output/$exonPositionsTranscriptPositions
 


#rm $exonPositionsSorted
python ./Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py\
 $output/Unique_DNA_Regions_for_riboseq.bed\
 $output/$exonPositionsTranscriptPositions\
 $output/Unique_DNA_Regions_genomic.bed

python ./Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py\
 $output/Unique_Protein_Regions_transcript_coords.bed\
 $output/$exonPositionsTranscriptPositions\
 $output/Unique_Protein_Regions_genomic.bed

# ----- Get the genomic positions for all valid Split-ORFs ----- #
python ./Genomic_scripts_18_10_24/create_ORF_coords_bed_file.py\
 $output/UniqueProteinORFPairs.txt\
  $output/ORF_transcript_coords.bed


python ./Genomic_scripts_18_10_24/genomic_DNA_regions_polars.py\
 $output/ORF_transcript_coords.bed\
 $output/$exonPositionsTranscriptPositions\
 $output/ORF_genomic_coords.bed


# ----- Test functionality of unique region to genomic coordinate conversion ----- #
python ./Genomic_scripts_18_10_24/test/test_transcriptomic_to_genomic_coordinates.py\
 ./Genomic_scripts_18_10_24/test/test_conversion_with_gen_trans.bed\
 $output/Unique_DNA_Regions_genomic.bed


# ----- calculate overlap between unique DNA and protein regions genomic         ----- #
echo "calculate genomic overlap between unique DNA and protein regions"
bedtools intersect\
 -a $output/Unique_DNA_Regions_genomic.bed\
 -b $output/Unique_Protein_Regions_genomic.bed\
 > $output/Unique_Regions_Overlap_genomic.bed


# ----- calculate overlap between unique DNA and protein regions transcriptomic         ----- #
echo "calculate transcriptomic overlap between unique DNA and protein regions"
bedtools intersect\
 -a $output/Unique_DNA_Regions_for_riboseq.bed\
 -b $output/Unique_Protein_Regions_transcript_coords.bed\
 > $output/Unique_Regions_Overlap_transcriptomic.bed



# ----- Create a report with basics statistics of the uniqueness scripts                       ----- #
R -e "rmarkdown::render('Extended_Pipeline_new.Rmd',output_file='$output/Uniqueness_Report.html',
params=list(args = c('/Output/run_$timestamp/Unique_DNA_Regions.fa', 
'/Output/run_$timestamp/Unique_Protein_Regions.fa',
'/Output/run_$timestamp/Unique_DNA_Regions_gt20_filtered.bed', 
'/Output/run_$timestamp/Unique_Protein_Regions_gt8_valid_filtered.bed',
'/Output/run_$timestamp/UniqueProteinORFPairs.txt',
'/Output/run_$timestamp/Unique_DNA_Regions_for_riboseq.bed',
'/Output/run_$timestamp/Unique_Protein_Regions_transcript_coords.bed',
'/Output/run_$timestamp/Unique_DNA_regions_genomic.bed',
'/Output/run_$timestamp/Unique_Protein_Regions_genomic.bed',
'/Output/run_$timestamp/Unique_Regions_Overlap_transcriptomic.bed',
'/Output/run_$timestamp/Unique_Regions_Overlap_genomic.bed')))"



# ----- get Split-ORF genes for GO analysis       ----- #
python ./SplitOrfs-master/get_so_genes_for_go_analysis.py \
 $output/UniqueProteinORFPairs.txt \
 $output/SplitOrfGeneFile.txt 




