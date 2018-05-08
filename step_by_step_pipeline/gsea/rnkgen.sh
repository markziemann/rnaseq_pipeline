#!/bin/bash
#rnkgen.sh converts a differential gene expression spreadsheet (XLS)
#into a rank file (RNK)

#Specify the input file
#XLS=$1
for XLS in *tsv ; do

#Specify the gene ID column
#ID=$2
ID=1

#Specify the fold change value column
#FC=$3
FC=2

#Specify the raw p-value column
#P=$4
P=5

#Specify ortholog maping
#ORTH=$5
ORTH=mouse2human.txt.sort

RNK=${XLS}.rnk
HUM=${RNK}.hum

sed 1d $XLS | tr -d '"' \
| awk -v I=$ID -v F=$FC -v P=$P '{FS="\t"} {print $I, $F, $P}' \
| awk '$2!="NA" && $3!="NA"' \
| awk '{s=1} $2<0{s=-1} {print $1"\t"s*-1*log($3)/log(10)}' \
| sort -k2gr \
| sed 's/inf$/305/' > $RNK

sed 's/_/\t/' $RNK \
| sort -k 1b,1 \
| join -1 2 -2 1 $ORTH - \
| awk '{OFS="\t"} {print $0,$5*$5}' \
| sort -k6gr \
| awk '!arr[$4]++' \
| awk '{OFS="\t"} !arr[$3]++ {print $3,$5}' \
| sort -k2gr > $HUM

sed -i -e '1iGeneID\tScore' $RNK
sed -i -e '1iGeneID\tScore' $HUM

done
