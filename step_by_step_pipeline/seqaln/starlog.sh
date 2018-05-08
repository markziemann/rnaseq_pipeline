#!/bin/bash

echo 'Number_of_input_reads_
Average_input_read_length_
Uniquely_mapped_reads_number_
Uniquely_mapped_reads_%_
Average_mapped_length_
Number_of_splices:_Total_
Number_of_splices:_Annotated_(sjdb)_
Number_of_splices:_GT/AG_
Number_of_splices:_GC/AG_
Number_of_splices:_AT/AC_
Number_of_splices:_Non-canonical_
Mismatch_rate_per_base,_%_
Deletion_rate_per_base_
Deletion_average_length_
Insertion_rate_per_base_
Insertion_average_length_
Number_of_reads_mapped_to_multiple_loci_
%_of_reads_mapped_to_multiple_loci_
Number_of_reads_mapped_to_too_many_loci_
%_of_reads_mapped_to_too_many_loci_
%_of_reads_unmapped:_too_many_mismatches_
%_of_reads_unmapped:_too_short_
%_of_reads_unmapped:_other_' > tmp.txt

for i in `cat tmp.txt`
do
echo -n "$i "
TAG=`echo $i | sed 's/_/ /g' `
grep -h "$TAG" *Log.final.out | cut -d '|' -f2 | sed 's/^[ \t]*//g;s/*[\t ]$//g' \
| tr '\n' '\t' | sed 's/$/\n/' | grep -v READS
done > Starlog.txt



for TAB in *ReadsPerGene.out.tab ; do
  NAME=$(echo $TAB | sed 's/.fq-trimmed.fastq_ReadsPerGene.out.tab//' )
  POS=$(tail -n+5 $TAB | cut -f3 | numsum)
  NEG=$(tail -n+5 $TAB | cut -f4 | numsum)
  RAT=$(echo $POS $NEG | awk '{print $1/$2}' | numround -n 0.001)
  ASS=$(tail -n+5 $TAB | cut -f4 | numsum)
  UNMAP=$(head -1 $TAB | cut -f4)
  MULTI=$(head -2 $TAB | tail -1 | cut -f4)
  NOFEAT=$(head -3 $TAB | tail -1 | cut -f4)
  AMBIG=$(head -4 $TAB | tail -1 | cut -f4)
  echo $NAME $RAT $ASS $UNMAP $MULTI $NOFEAT $AMBIG
done

