library("plyr")
library("statmod")
library("edgeR")
library("locfit")
library("reshape2")
library("parallel")
library("gplots")

tmp<-read.table("3col.tsv",header=F)

x<-as.matrix(acast(tmp, V2~V1, value.var="V3"))
colnames(x)<-sub(".fq-trimmed.fastq", "", colnames(x))

#dont forget gene names
g<-read.table("Mus_musculus.GRCm38.91.gnames.txt",row.names=1)
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

colfunc <- colorRampPalette(c("blue", "white", "red"))

#################################################
# Do some analysis of all samples MDS plot
##################################################
y<-x[which(rowSums(x)/ncol(x)>=(10)),]
z<-log2(sweep(y,2,colSums(y),'/')*1000000)
z <- z[!is.infinite(rowSums(z)),]

pdf("MDSplots.pdf")
plot(cmdscale(dist(t(z))), xlab="Coordinate 1", ylab="Coordinate 2", type = "n",main="MDS plot all samples")
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
y <- estimateDisp(y, design,robust=TRUE,prior.df=20)
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
rnk<-as.data.frame(-log10(dge$PValue+1E-300)/sign(dge$logFC))
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
#heatmap(as.matrix(z),margins = c(6,16))
heatmap.2(as.matrix(z),col=colfunc(25),scale="row", trace="none",margins = c(8,16), cexRow=.5, main="top significant genes" ,  keysize=1)
dev.off()

# gene barchart
rownames(dge)=dge$Row.names
pdf(file=paste(PFX,"_genes.pdf",sep=""))
for (i in 1:20) {
COLS=grep(".x",colnames(dge))
b<-t(dge[i,COLS])
logFC=dge[i,grep("logFC",colnames(dge)) ]
Pval=dge[i,grep("PVal",colnames(dge)) ]
FDR=dge[i,grep("FDR",colnames(dge)) ]
HEADER=colnames(b)
SUBHEAD=paste("logFC=",signif(logFC,3)," Pval=",signif(Pval,3)," FDR=",signif(FDR,3),sep="")
par(mar=c(10,8,4,5))
barplot(b,beside=T,main=HEADER,names.arg=sub(".x","",rownames(b)),las=2,ylab="CPM",cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1)
mtext(SUBHEAD)
}
dev.off()


}

#mapply(dge_cycle,as.vector(a$Var1),as.vector(a$Var2))
mcmapply(dge_cycle,as.vector(a$Var1),as.vector(a$Var2),mc.cores=32)

