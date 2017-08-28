#!/bin/bash
THISPATH=$(pwd)
#RNAseq pipeline from start to finish

###################################################
# User defined parameters.
# These need to be configured for each analysis
###################################################

# Specify project name
PROJ_PATH=$THISPATH/cyp_pfi

# Specify user name and password
USER=JohnDoe
PW=secret

# Specify reference genome and gtf
GENOME_URL=ftp://ftp.ensembl.org/pub/release-90/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
GTF_URL=ftp://ftp.ensembl.org/pub/release-90/gtf/homo_sapiens/Homo_sapiens.GRCh38.90.gtf.gz

# Specify strandedness of the library prep
# 0=unstranded, 1=pos strand, 2=neg strand
STRAND=2

# Specify path to folder containing gsea JAR file GMT files
GSEA_UTILS=/path/to/gsea/utils/

# Specify URL for data download
# must be a directory
URL=ftp://ngs_service/NGS_RUN_ABC123


###################################################
# Don't modify the following text. it should just work.
###################################################

REF_DIR=$PROJ_PATH/refgenome
FASTQ_DIR=$PROJ_PATH/fastq
ALN_DIR=$PROJ_PATH/aln
EDGER_DIR=$PROJ_PATH/edger
GSEA_DIR=$PROJ_PATH/gsea
INTEGRATION_DIR=$PROJ_PATH/integration

mkdir -p $PROJ_PATH $REF_DIR $FASTQ_DIR $ALN_DIR $EDGER_DIR $GSEA_DIR $INTEGRATION_DIR
cd $PROJ_NAME


###################################################
# download fastq
###################################################
cd $FASTQ_DIR

wget -r --no-host-directories --ftp-user=$USER --ftp-password=$PW "$URL"
ftp://melbourne-ftp.agrf.org.au/AGRF_CAGRF15810_CB6N8ANXX
RAWDATA_DIR=$FASTQ_DIR/$(ls)
echo $RAWDATA_DIR
cd $RAWDATA_DIR
#cat checksums.md5 | parallel --pipe -N1 md5sum -c | tee checksum_report.txt

###################################################
# download and index reference genomes
###################################################
cd $REF_DIR

if [ ! -r "SAindex" ] ; then
  echo reference not present. fetching now
  axel -n 5 $GENOME_URL
  axel -n 5 $GTF_URL
  pigz -d *gz

fi

  GNAMES=$(echo $GTF | sed 's#.gtf#.gnames.txt#')
  GTF=$(find . | grep gtf$)
  FA=$(find . | egrep '(fa$|.fna$)')
  CWD=$(pwd)

if [ ! -r "SAindex" ] ; then

  STAR --runMode genomeGenerate \
  --sjdbGTFfile $GTF \
  --genomeDir $CWD  \
  --genomeFastaFiles $CWD/$FA \
  --runThreadN $(nproc)

  grep -w gene $GTF | cut -d '"' -f2,6 | tr '"' '\t' | sort -k 1b,1 > $GNAMES

fi


###################################################
# Trim reads with skewer
###################################################
cd $ALN_DIR
#ln $RAWDATA_DIR/*fastq.gz .
rm e_*fastq.gz
trim(){
FQZ=$1
skewer -l 20 -q 20 -t $(nproc) $FQZ
}
export -f trim
#parallel -j8 trim ::: *fastq.gz

###################################################
# Map reads to the genome with star
###################################################
#STAR --genomeLoad LoadAndExit --genomeDir $REF_DIR
runstar(){
FQ=$1
REF=$2
STAR --runThreadN $(nproc) --quantMode GeneCounts --genomeLoad LoadAndKeep \
--outFileNamePrefix ${FQ}_ --outSAMtype None --genomeDir $REF --readFilesIn=$FQ
}
export -f runstar
#parallel -j 3 runstar ::: *fastq ::: $REF_DIR
#STAR --genomeLoad Remove --genomeDir $REF_DIR


###################################################
# Extract counts into a 3 colum table sample-gene-count
###################################################
>3col.txt

CUT_COL=$((STRAND+2))
for CNT in *ReadsPerGene.out.tab ; do
  NAME=$(echo $CNT | cut -d '_' -f-2)
  tail -n+5 $CNT | cut -f1,$CUT_COL | sed "s/^/${NAME}\t/" >> 3col.txt
done

cp 3col.txt $REF_DIR/$GNAMES $EDGER_DIR

cd $EDGER_DIR
cat <<'EOF' > edgeR_script.R
library("plyr")
library("statmod")
library("edgeR")
library("locfit")
library("reshape2")
library("parallel")

tmp<-read.table("3col.txt",header=F)

x<-as.matrix(acast(tmp, V2~V1, value.var="V3"))
#dont forget gene names
g<-read.table("Homo_sapiens.GRCh38.90.gnames.txt",row.names=1)
x<-merge(g,x,by=0)
rownames(x)=paste(x$Row.names,x$V2,sep="_")
x$Row.names=NULL
x$V2=NULL

samplesheet<-as.data.frame(colnames(x))
colnames(samplesheet)="sample"
samplesheet$grp<-sub("_.*","",samplesheet$sample)
rownames(samplesheet)=samplesheet$sample
samplesheet$sample=NULL

a<-expand.grid(unique(samplesheet$grp),unique(samplesheet$grp))
a<-subset(a,a$Var1!=a$Var2)
a<-a[order(a$Var1,a$Var2),]

dgelist=NULL

#################################################
# Do some analysis of all samples MDS plot
##################################################
y<-x[which(rowSums(x)/ncol(x)>=(10)),]
z<-log2(sweep(y,2,colSums(y),'/')*1000000)
z <- z[!is.infinite(rowSums(z)),]

pdf("MDSplot.pdf")
plot(cmdscale(dist(t(z))), xlab="Coordinate 1", ylab="Coordinate 2", type = "n")
text(cmdscale(dist(t(z))), labels=colnames(z), )
dev.off()

#################################################
#define the dge pipeline
#################################################
dge_cycle<-function(GROUP1,GROUP2){
message(paste(GROUP1,GROUP2))

PFX=paste(GROUP1,"vs",GROUP2,sep="_")
# format design matrix for DGE without pairing information
des<-subset(samplesheet,grp==GROUP1 | grp==GROUP2)
group<-factor(des$grp, levels=c(GROUP1, GROUP2), ordered=TRUE)
design<-model.matrix(~group)
rownames(design)=rownames(des)

#subset counts
x1<-as.data.frame(x[,rownames(design)])
x1<-x1[which(rowSums(x1)/ncol(x1)>=(10)),]

#MDS PLOT
pdf(paste(PFX,"_MDSplot.pdf",sep=""))
plot(cmdscale(dist(t(x1))), xlab="Coordinate 1", ylab="Coordinate 2", type = "n")
text(cmdscale(dist(t(x1))), labels=colnames(x1), )
dev.off()

# run the dge analysis
y<-DGEList(counts=x1)
y <- calcNormFactors(y)
y <- estimateDisp(y, design,robust=TRUE,prior.df=1)
fit <- glmFit(y, design)
lrt <- glmLRT(fit)
dge<-as.data.frame(topTags(lrt,n=1000000))
dge$dispersion<-lrt$dispersion
dge<-merge(dge,lrt$fitted.values,by='row.names')
rownames(dge)=dge$Row.names
dge$Row.names=NULL
dge<-merge(dge,y$counts,by='row.names')
dge<-dge[order(dge$PValue),]
write.table(dge,file=paste(PFX,".tsv",sep=""),sep="\t",quote=F,row.names=F)
dgelist<-rbind(dgelist,as.data.frame(subset(dge,FDR<0.05,select=c("Row.names"))))
head(dge)

##output rank file
rnk<-as.data.frame(sign(dge$logFC)/log10(dge$PValue+1E-300))
colnames(rnk)="Score"
rownames(rnk)=dge$Row.names
write.table(rnk,file=paste(PFX,".rnk",sep=""),sep="\t",quote=F,row.names=T)

# MA plot
dge2<-subset(dge,FDR<0.05)
cntdge=nrow(dge2)
cntup=nrow(subset(dge2,logFC>0))
cntdn=nrow(subset(dge2,logFC<0))
HEADING=paste(PFX, cntdge, "DGEs", cntup, "up", cntdn, "dn")
pdf(file=paste(PFX,"_MAplot.pdf",sep=""))
plot(dge$logCPM,dge$logFC,main=HEADING,col="gray",pch=19,cex=0.5,xlab="log2(CPM)",ylab="log2(Fold Change)")
points(dge2$logCPM,dge2$logFC,col="red",pch=19,cex=0.5)
#text(y$logCPM+1,y$logFC,labels=rownames(y),col="black",pch=19,cex=0.6)
dev.off()

#TODO heatmap
cols=paste(grep(".x", names(dge), value = TRUE))
z<-as.data.frame(dge[,cols])
rownames(z)=dge$Row.names
z<-z[1:50,]
pdf(paste(PFX,"_heatmap.pdf",sep=""),width=7, height=8)
heatmap(as.matrix(z))
dev.off()

}

#mapply(dge_cycle,as.vector(a$Var1),as.vector(a$Var2))
mcmapply(dge_cycle,as.vector(a$Var1),as.vector(a$Var2),mc.cores=32)
EOF
Rscript edgeR_script.R


#############################################################
echo start GSEA
#############################################################

cp *rnk $GSEA_DIR
cd $GSEA_DIR
sed -i 's/Score/Accession_GeneID\tScore/' *rnk
for RNK in *rnk ; do
  sed 1d $RNK | cut -d '_' -f2- | awk '{OFS="\t"} {print $1,$2,$2*$2}' \
  | sort -k3gr | awk '{OFS="\t"} !arr[$1]++ {print $1,$2}' > tmp
  sed -e '1i\GeneID\tScore' tmp > $RNK
done

cp $GSEA_UTILS/* .
run_gsea(){
GSEAJAR=*jar
RNK=$1
GMT=$2
echo $GSEAJAR $RNK $GMT
java -Xmx4096m -cp $GSEAJAR xtools.gsea.GseaPreranked  \
-gmx $GMT -collapse false -mode Max_probe \
-norm meandiv -nperm 1000 -rnk $RNK -scoring_scheme classic \
-rpt_label ${RNK}.${GMT} -include_only_symbols true -make_sets true \
-plot_top_x 20 -rnd_seed timestamp -set_max 10000 -set_min 10 -zip_report false \
-out . -gui false
}
export -f run_gsea
parallel -j8 run_gsea ::: *.rnk ::: *.gmt

echo 'GeneSetName	GeneSetSize	ES	NES	p-val	FDR	FWER' > header.txt

for GSEADIR in `ls | grep GseaPreranked | grep -v xls$` ; do
  awk '{FS="\t"} {OFS="\t"} $8<0.05 {print $1,$4,$5,$6,$7,$8,$9} ' $GSEADIR/gsea_report_for_na_*xls \
  | cat header.txt - > $GSEADIR.xls
done


cat <<'EOF' >> plotGSEA.R
path = "."
infiles <- dir(pattern='\\.xls$')
plot.GSEA <- function(file){
  x <-read.table(file,header=T,row.names=1)
  y<-head(x[order(x$ES),],n=20L)
  y<-y[which(y$ES<0),]
  z<-head(x[order(-x$ES),],n=20L)
  z<-z[which(z$ES>0),]
  df <- rbind(y,z)
  df<-df[order(df$ES),]
  barplot(df$ES,main=file,xlab="GSEA enrichment score",horiz=TRUE,names.arg=row.names(df))
}
pdf(file="gsea_results.pdf",width=15,height=10)
par(las=2) ; par(mar=c(10,50,1,1))
lapply(infiles,plot.GSEA)
dev.off()
lapply(infiles,plot.GSEA)
EOF
Rscript plotGSEA.R









