#!/bin/bash
echo 'GeneSetName	GeneSetSize	ES	NES	p-val	FDR	FWER' > header.txt

for GSEADIR in `ls | grep GseaPreranked | grep -v xls$` ; do
  awk '{FS="\t"} {OFS="\t"} $8<0.05 {print $1,$4,$5,$6,$7,$8,$9} ' $GSEADIR/gsea_report_for_na_*xls \
  | cat header.txt - > $GSEADIR.xls
done


