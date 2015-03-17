#!/bin/bash

bamfile1=/home/rode/sead/projects/PilotJan2015/mapping_PilotJan2015/mapping_PilotJan2015_chr1.sorted.bam

bamfile2=/home/rode/sead/projects/PilotJan2015/realign_PilotJan2015/mapping_PilotJan2015_realignedchr1.sorted.bam

bamfilename1=$(basename $bamfile1 .bam)
bamfilename2=$(basename $bamfile2 .bam)


if [[ -e ${bamfilename1}.sorted.bam ]];then

	samtools sort -n $bamfilename1 $bamfilename1.sorted

fi

if [[ -e ${bamfilename2}.sorted.bam ]];then

	samtools sort -n $bamfilename2 $bamfilename2.sorted

fi


if [[ -e tmp5 ]];then

	rm tmp5
        echo "Former folder (tmp5) has been removed."

fi




#samtools view $bamfile1 | grep -o "NM:i:.." | awk -F":" 'BEGIN {s=0}{s+=$3}END{print s}'

samtools view $bamfilename1.sorted.bam |head -n 10000 | cut -f6  > tmp1
samtools view $bamfilename2.sorted.bam |head -n 10000 | cut -f6  > tmp2

paste tmp1 tmp2 > tmp3

awk '{if ($1!=$2) {print NR} }' tmp3 > tmp4


for line in $(cat tmp4);do

	samtools view $bamfilename1.sorted.bam |head -n 10000 | sed -n ${line}p >> tmp5
	samtools view $bamfilename2.sorted.bam |head -n 10000 | sed -n ${line}p >> tmp5

done


