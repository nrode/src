#!/bin/bash

#Usage: Extract the name of each R1 sequence in a fastq.gz file and compute the length of the corresponding R1 and R2 sequences
#input: name of the fastq.gz file with R1 sequences

fastqfile1=$1
fastqfile2=${fastqfile1/R1/R2}
namefile=$(basename $fastqfile1 .fastq.gz)
numseq=$(($(zcat $fastqfile1 | grep -c .)/4+1))
#numseq=3

for ((i=1;i<$numseq;i=i+1));do

	namerow=$(((i-1)*4+1))
	seqrow=$(((i-1)*4+2))

	nameseq=$(zcat $fastqfile1 | sed -n ${namerow}p )
	len1=$(zcat $fastqfile1 | sed -n ${seqrow}p | awk '{print length;}')
	len2=$(zcat $fastqfile2 | sed -n ${seqrow}p | awk '{print length;}')

	echo $nameseq $len1 $len2 $(($len1-$len2)) | gzip --stdout >> $namefile.count.gz

done
