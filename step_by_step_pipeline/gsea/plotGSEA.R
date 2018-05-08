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
