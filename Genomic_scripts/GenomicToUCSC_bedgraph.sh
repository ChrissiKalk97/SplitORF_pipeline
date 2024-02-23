in=$1
ucOut=$(echo $in | rev | cut -f 2- -d '.' | rev)_UCSC.bed
scafOut=$(echo $ucOut | rev | cut -f 2- -d '.' | rev)_no_scaffold.bed

source ./Genomic_scripts/ScaffoldUmrechnung.sh $in
./Genomic_scripts/chromToUcsc -a ./Genomic_scripts/chromAlias_38p14_18_12_23.tsv -i $in -o $ucOut
cat $ucOut | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3,$4; else print $1,$2,$3,$4}' > $scafOut