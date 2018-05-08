#!/bin/bash
GTF=$(find . | grep gtf$)
FA=$(find . | egrep '(fa$|.fna$)')
CWD=$(pwd)
GNAMES=$(echo $GTF | sed 's#.gtf#.gnames.txt#')

  STAR --runMode genomeGenerate \
  --sjdbGTFfile $GTF \
  --genomeDir $CWD  \
  --genomeFastaFiles $CWD/$FA \
  --runThreadN $(nproc)

  grep -w gene $GTF | cut -d '"' -f2,6 | tr '"' '\t' | sort -k 1b,1 > $GNAMES

