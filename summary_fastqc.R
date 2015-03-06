#!/bin/Rscript
args<-commandArgs(trailingOnly=TRUE)
dir1=args[1]
dir2=args[2]

print(dir2)
setwd(paste("~/sead/projects/",dir1,"/",dir2,sep=''))
source("/home/rode/src/quantgen/utils_fastqc.R")

results<-read.fastq.zips()
nreads<-nreads.fastqc(results)
write.table(nreads,"readcounts.txt",quote=F,col.names=F)

pdf(file=paste("FastqcSummary_",dir2,".pdf",sep=""))
#Barplot with the number of reads per fastq.gz file
barplot.nreads.fastqc(nreads)
#Histogram with the distribution of the phred quality
distribquality<-quals.fastqc(results,perc=T,nreads=nreads)
plot.nbseq.qual(distribquality,perc=T,main="Phred distribution",legend.cex=0.7)
#Plot the percentage of adapter content across base positions
nbadp <- names(results[[1]][["Adapter_Content"]])[-1]
for (i in nbadp){
percadaptercontent<-adp.contents.fastqc(results,adp=i)
namegraph<-paste("Adapter content",i,sep=" ")
plot.adp.content(percadaptercontent,max.dataset.per.plot=25,main=namegraph,legend.cex=0.7)
}
#Plot with seq length numbers
distriblength<-seq.lengths.fastqc(results)
distriblength<-log10(distriblength+1)
plot.seq.length(distriblength,main="Sequence Length Distribution All",ylab="log10(Number of sequences)",legend.cex=0.7)
distriblength<-seq.lengths.fastqc(results)
distriblength<-log10(distriblength+1)
plot.seq.length(distriblength,max.dataset.per.plot=2, main="Sequence Length Distribution",ylab="log10(Number of sequences)",legend.cex=0.7)
write.table(distriblength,"seq.length.distrib.txt",quote=F)
dev.off()

