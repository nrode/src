#!/bin/bash


#Split a bam file according to the different read groups included in the header


inbam=$1
inbamname=$(basename $1 .bam)

#index=$2

index="RPI2|RPI3|RPI5"
#index="RPI6|RPI7|RPI9"

loc="scaffold0002"
out=${3:-bamlist}

if [[ -e $out ]];then

	rm $out
	echo "Former list ("${out}") has been removed"

fi

RG=$(samtools view -H $inbam | grep "^@RG" |cut -f2|cut -f2 -d':'| grep -v "Mix" | egrep $index)

set -xv

for ind in $RG;do

	filename=${inbamname}_${ind}.bam
	samtools view -b -r $ind $inbam $loc > $filename
#	samtools view -H -b -r $ind $inbam $loc > $filename
	echo $filename >> $out

done

