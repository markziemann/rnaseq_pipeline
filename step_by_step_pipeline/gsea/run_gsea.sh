#!/bin/bash
run_gsea(){
RNK=$1
GMT=$2
#echo /data/projects/mziemann/app/gsea2-2.2.2.jar $RNK $GMT
#gsea2-3.0_beta_2.jar
GSEAJAR=gsea-3.0_beta_3.jar
echo $GSEAJAR $RNK $GMT
java -Xmx4096m -cp $GSEAJAR xtools.gsea.GseaPreranked  \
-gmx $GMT -collapse false -mode Max_probe \
-norm meandiv -nperm 1000 -rnk $RNK -scoring_scheme classic \
-rpt_label ${RNK}.${GMT} -include_only_symbols true -make_sets true \
-plot_top_x 20 -rnd_seed timestamp -set_max 10000 -set_min 10 -zip_report false \
-out . -gui false
}
export -f run_gsea
parallel run_gsea ::: *.hum ::: *.gmt

exit
GSEA
java -Xmx512m xtools.gsea.GseaPreranked -gmx c2.cp.reactome.v5.0.symbols.gmt \
-norm meandiv -nperm 1000 -rnk pairededgeR_xlshum.rnk -scoring_scheme classic \
-rpt_label gsea3_for_SVG -create_gcts false -create_svgs true -help false \
-make_sets true -plot_top_x 50 -rnd_seed timestamp -set_max 5000 -set_min 5 \
-zip_report false -out gsea -gui false

