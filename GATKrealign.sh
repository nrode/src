#!/bin/bash

project_folder=$1
loc="chr1"
cd ~/sead/projects/${project_folder}

realdir="realign_${project_folder}"

if [[ -d $realdir ]]; then
        echo "Former folder ($realdir) has been removed."
#        rm -r $realdir;
fi

#mkdir $realdir

cd $realdir


#Look for positions to be realigned
#Without encoding option
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -o ./RTC.intervals -I ../mapping_${project_folder}/mapping_${project_folder}.bam -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta


#Perform the realignment
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ../mapping_${project_folder}/mapping_${project_folder}.bam -L $loc -targetIntervals ./RTC.intervals -o ./mapping_${project_folder}_realigned${loc}.bam

#Index realignement
#samtools index realign_PilotJan2015/mapping_PilotJan2015_realigned.bam


