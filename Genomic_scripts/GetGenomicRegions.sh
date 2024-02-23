#!/bin/bash

#Help message:
usage="
Usage: ./GetGenomicRegions.sh [-h] GenomicAnnotation.bed UniqueDNA.bed UniqueDNAOut_with_chr.bed UniqueDNAOut_with_trID.bed UniqueProtein.bed UniqueProteinOut_with_chr.bed UniqueProteinOut_with_trID.bed

GenomicAnnotation.bed 			a BED file containing the genomic positions for all transcripts
ExonAnnotation.bed			a Bed file containing the genomic positions for all exons
UniqueDNA.bed 				the unique_DNA_regions_merged.bed file from the main pipeline
UniqueProtein.bed			the unique_protein_regions_merged_valid.bed file from the main pipeline
ORFProtein.bed			the Valid_ORF_Proteins.bed file from the main pipeline

where:
-h	show this help"

#available options for the programm
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

#Check that all arguments are provided and are in the correct file format. Give error messages according to errors in the call and exit the programm if something goes wrong
RED='\033[0;31m' #Red colour for error messages
NC='\033[0m'	 #No colour for normal messages

#check for right number of arguments
if [ "$#" -ne 5 ]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
#check if any argument is a directory
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ] || [ -d "$5" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
#check if every argument has the correct file extension
elif ! [[ $1 =~ \.bed$ ]] || ! [[ $2 =~ \.bed$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.bed$ ]] || ! [[ $5 =~ \.bed$ ]]; then
  echo -e "${RED}
ERROR while executing!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi
GenomicAnnotation=$1
ExonAnnotation=$2
UniqueDNA=$3
UniqueDNAcorrect=$(echo $UniqueDNA | rev | cut -f 2- -d '.' | rev)_unspliced_positions.bed
UniqueDNAOutChr=$(echo $UniqueDNA | rev | cut -f 2- -d '.' | rev)_genomic_with_chr.bed
UniqueDNAOutTrID=$(echo $UniqueDNA | rev | cut -f 2- -d '.' | rev)_genomic_with_trID.bed
UniqueProtein=$4
UniqueProteincorrect=$(echo $UniqueProtein | rev | cut -f 2- -d '.' | rev)_unspliced_positions.bed
UniqueProteinOutChr=$(echo $UniqueProtein| rev | cut -f 2- -d '.' | rev)_genomic_with_chr.bed
UniqueProteinOutTrID=$(echo $UniqueProtein | rev | cut -f 2- -d '.' | rev)_genomic_with_trID.bed
OrfProteins=$5
OrfProteinscorrect=$(echo $OrfProteins | rev | cut -f 2- -d '.' | rev)_unspliced_positions.bed
OrfOutChr=$(echo $OrfProteins | rev | cut -f 2- -d '.' | rev)_genomic_with_chr.bed
OrfOutTrID=$(echo $OrfProteins | rev | cut -f 2- -d '.' | rev)_genomic_with_trID.bed

python ./Genomic_scripts/AddExonPositions.py $ExonAnnotation $UniqueDNA $UniqueDNAcorrect
python ./Genomic_scripts/genomicpositions.py $GenomicAnnotation $UniqueDNAcorrect $UniqueDNAOutChr $UniqueDNAOutTrID
source ./Genomic_scripts/GenomicToUCSC.sh $UniqueDNAOutChr

python ./Genomic_scripts/AddExonPositions_protein.py $ExonAnnotation $UniqueProtein $UniqueProteincorrect
python ./Genomic_scripts/genomicpositions_protein.py $GenomicAnnotation $UniqueProteincorrect $UniqueProteinOutChr $UniqueProteinOutTrID
source ./Genomic_scripts/GenomicToUCSC.sh $UniqueProteinOutChr

python ./Genomic_scripts/AddExonPositions_ORF.py $ExonAnnotation $OrfProteins $OrfProteinscorrect
python ./Genomic_scripts/genomicOrfPos.py $GenomicAnnotation $OrfProteinscorrect $OrfOutChr $OrfOutTrID
source ./Genomic_scripts/GenomicToUCSC.sh $OrfOutChr
