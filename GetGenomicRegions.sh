#!/bin/bash

#Help message:
usage="
Usage: ./GetGenomicRegions.sh [-h] GenomicAnnotation.bed UniqueDNA.bed UniqueDNAOut_with_chr.bed UniqueDNAOut_with_trID.bed UniqueProtein.bed UniqueProteinOut_with_chr.bed UniqueProteinOut_with_trID.bed

GenomicAnnotation.bed 			a BED file containing the genomic positions for all used transcripts
UniqueDNA.bed 				the unique_DNA_regions_merged.bed file from the main pipeline
UniqueDNAOut_with_chr.bed 			the name for the unique DNA regions file with genomic positions (with chromosome/scaffold name)
UniqueDNAOut_with_trID.bed 			the name for the unique DNA regions file with genomic positions (with geneID|transcriptID)
UniqueProtein.bed			the unique_protein_regions_merged_valid.bed file from the main pipeline
UniqueProteinOut_with_chr.bed 			the name for the unique protein regions file with genomic positions (with chromosome/scaffold name)
UniqueProteinOut_with_trID.bed 			the name for the unique protein regions file with genomic positions (with geneID|transcriptID)

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
if [ "$#" -ne 10 ]; then #check for right number of arguments
  echo -e "${RED}
ERROR while executing the Pipeline!
Wrong number of arguments.${NC}"
  echo "$usage" >&2
  exit 1
#check if any argument is a directory
elif  [ -d "$1" ] || [ -d "$2" ] || [ -d "$3" ] || [ -d "$4" ] || [ -d "$5" ] || [ -d "$6" ] || [ -d "$7" ] || [ -d "$8" ] || [ -d "$9" ] || [ -d "$10" ]; then
  echo -e "${RED}
ERROR while executing the Pipeline!
One or more of the arguments are directories.${NC}"
  echo "$usage" >&2
  exit 1
#check if every argument has the correct file extension
elif ! [[ $1 =~ \.bed$ ]] || ! [[ $2 =~ \.bed$ ]] || ! [[ $3 =~ \.bed$ ]] || ! [[ $4 =~ \.bed$ ]] || ! [[ $5 =~ \.bed$ ]] || ! [[ $6 =~ \.bed$ ]] || ! [[ $7 =~ \.bed$ ]]|| ! [[ $8 =~ \.bed$ ]] || ! [[ $9 =~ \.bed$ ]] || ! [[ $10 =~ \.bed$ ]]; then
  echo -e "${RED}
ERROR while executing!
One or more of the arguments are not in the specified file format.${NC}"
  echo "$usage" >&2
  exit 1
fi
GenomicAnnotation=$1
UniqueDNA=$2
UniqueDNAOutChr=$3
UniqueDNAOutTrID=$4
UniqueProtein=$5
UniqueProteinOutChr=$6
UniqueProteinOutTrID=$7
OrfProteins=$8
OrfOutChr=$9
OrfOutTrID=$10

python ./Uniqueness_scripts/genomicpositions.py $GenomicAnnotation $UniqueDNA $UniqueDNAOutChr $UniqueDNAOutTrID
source ./ScaffoldUmrechnung.sh $UniqueDNAOutChr
UniqueDNAUCSC=$(echo $UniqueDNAOutChr | rev | cut -f 2- -d '.' | rev)_UCSC.bed
./chromToUcsc -a ./hg38.chromAlias.tsv -i $UniqueDNAOutChr -o $UniqueDNAUCSC
UniqueDNAnoScaffolds=$(echo $UniqueDNAUCSC | rev | cut -f 2- -d '.' | rev)_no_scaffolds.bed
cat $UniqueDNAUCSC | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3; else print $1,$2,$3}' > $UniqueDNAnoScaffolds

python ./Uniqueness_scripts/genomicpositions_protein.py $GenomicAnnotation $UniqueProtein $UniqueProteinOutChr $UniqueProteinOutTrID
source ./ScaffoldUmrechnung.sh $UniqueProteinOutChr
UniqueProteinUCSC=$(echo $UniqueProteinOutChr | rev | cut -f 2- -d '.' | rev)_UCSC.bed
../../chromToUcsc -a ../../hg38.chromAlias.tsv -i $UniqueProteinOutChr -o $UniqueProteinUCSC
UniqueProteinnoScaffolds=$(echo $UniqueProteinUCSC | rev | cut -f 2- -d '.' | rev)_no_scaffolds.bed
cat $UniqueProteinUCSC | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3; else print $1,$2,$3}' > $UniqueProteinnoScaffolds

python ./Uniqueness_scripts/genomicOrfPos.py $GenomicAnnotation $OrfProteins $OrfOutChr $OrfOutTrID
source ./ScaffoldUmrechnung.sh $OrfOutChr
OrfUCSC=$(echo $OrfOutChr | rev | cut -f 2- -d '.' | rev)_UCSC.bed
../../chromToUcsc -a ../../hg38.chromAlias.tsv -i $OrfOutChr -o $OrfUCSC
OrfnoScaffolds=$(echo $OrfUCSC | rev | cut -f 2- -d '.' | rev)_no_scaffolds.bed
cat $OrfUCSC | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3; else print $1,$2,$3}' > $OrfnoScaffolds
