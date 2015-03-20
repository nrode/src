#!/bin/bash
# usage: ~/rode/sead/src/loop_bamsplitRG.sh project_folder  
# arguments: $1 = name of the project's folder

project_folder=$1
cd ~/sead/projects/$project_folder

dir=bamsplit_${project_folder}

if [[ -d $dir ]]; then
	echo "Former folder ("$dir") has been removed."
	rm -r $dir
fi

mkdir $dir

cd $dir
set -xv
indexlist="RPI2|RPI3 RPI6|RPI7 RPI4 RPI5 RPI8 RPI9"
#indexlist="RPI4"
for index in $indexlist;do

	../../../src/bamsplitRG.sh -b ../realign_${project_folder}/mapping_PilotJan2015_filtered_realigned.bam -i $index

done

indexA17list="RPI5 RPI9"
indexA17list="RPI5"
s=4
for indexA17 in $indexA17list; do
	samtools view -H mapping_PilotJan2015_filtered_realigned_$indexA17.bam | grep "^@RG" |cut -f2|cut -f2 -d':'| grep -v "A17" > indexlist.tmp

	samtools view -h -b -R indexlist.tmp mapping_PilotJan2015_filtered_realigned_$indexA17.bam > mapping_PilotJan2015_filtered_realigned_${indexA17}_Mix.bam	
	samtools view -h -b -r ${indexA17}_S${s}_A17 mapping_PilotJan2015_filtered_realigned_$indexA17.bam > mapping_PilotJan2015_filtered_realigned_${indexA17}_A17.bam
	s+=1
	rm indexlist.tmp 
	rm mapping_PilotJan2015_filtered_realigned_$indexA17.bam 
done

#samtools index *.bam
