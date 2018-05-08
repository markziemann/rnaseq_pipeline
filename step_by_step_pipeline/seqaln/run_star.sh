#!/bin/bash
REF=../refgenome

run_pipeline(){
CWD=$(pwd)
FQZ1=$1
REF=$2
FQZ2=$(echo $FQZ1 | sed 's/_1/_2/')
FQ1=$(echo $FQZ1 | sed 's/.gz$/-trimmed-pair1.fastq/')
FQ2=$(echo $FQZ1 | sed 's/.gz$/-trimmed-pair2.fastq/')
BASE=$(echo $1 | sed 's/.fq.gz//')
BAM=$BASE.bam

skewer -t $(nproc) -q 20 $FQZ1 $FQZ2

STAR --runThreadN $(nproc) --quantMode GeneCounts --genomeLoad LoadAndKeep  \
 --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000  \
 --genomeDir $REF --readFilesIn=$FQ1 $FQ2 --outFileNamePrefix $BASE.

rm $FQ1 $FQ2
}
export -f run_pipeline

parallel -j1 run_pipeline ::: *_1*.fq.gz ::: $REF
STAR --genomeLoad Remove --genomeDir $REF

for TAB in *ReadsPerGene.out.tab ; do
  tail -n +5 $TAB | cut -f1,4 | sed "s/^/${TAB}\t/"
done > 3col.tsv

