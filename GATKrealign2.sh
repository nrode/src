#!/bin/bash

project_folder=$1


#Look for positions to be realigned
#Without encoding option
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -o realign_PilotJan2015/RTC.intervals -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta


#Perform the realignment
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -targetIntervals realign_PilotJan2015/real/RTC.intervals -o realign_PilotJan2015/mapping_PilotJan2015_realigned.bam

#Index realignement
samtools index realign_PilotJan2015/mapping_PilotJan2015_realigned.bam


