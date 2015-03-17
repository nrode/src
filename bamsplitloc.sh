#!/bin/bash


#Split a bam file according to the user-defined position

inbam=$1
loc=$2

filename=$(basename $inbam .bam)_${loc}.bam

if [[ -e $out ]];then

	rm $out
	echo "Former list ("${out}") has been removed"

fi

	samtools view -b $inbam $loc > $filename


