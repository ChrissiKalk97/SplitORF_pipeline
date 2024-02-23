in=$1
ucOut=$(echo $in | rev | cut -f 2- -d '.' | rev)_UCSC.bed
scafOut=$(echo $ucOut | rev | cut -f 2- -d '.' | rev)_no_scaffold.bed

source ./Genomic_scripts/ScaffoldUmrechnung.sh $in
./Genomic_scripts/chromToUcsc -a ./hg38.chromAlias.tsv -i $in -o $ucOut
cat $ucOut | awk '{if(index($1,"_")>0)print substr($1,1,index($1,"_")-1),$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12; else print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > $scafOut