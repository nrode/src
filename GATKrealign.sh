#!/bin/bash

project_folder=$1
#loc="chr1"

cd ~/sead/projects/${project_folder}

realdir="realign_${project_folder}chr1"

if [[ -d $realdir ]]; then
        echo "Former folder ($realdir) has been removed."
        rm -r $realdir;
fi

mkdir $realdir

cd $realdir
pathtoinbam=../mapping_${project_folder}
inbam=mapping_${project_folder}_chr1.bam
outfiltered=$(basename ${inbam} .bam)_filtered.bam
outrealigned=$(basename ${outfiltered} .bam)_realigned.bam

#Filter on quality
samtools view -b -h -q 30 ${pathtoinbam}/${inbam} > $outfiltered

#Index realignement
samtools index ./${outfiltered}

#Look for positions to be realigned
#Without encoding option
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -o ./RTC.intervals -I ./${outfiltered} -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta

ls

#Perform the realignment
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ./${outfiltered} -targetIntervals ./RTC.intervals -o ./${outrealigned}

#Index realignement
#samtools index ./${outrealigned}



#Realignement for a given locus
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ../mapping_${project_folder}/mapping_${project_folder}.bam -L $loc -targetIntervals ./RTC.intervals -o ./mapping_${project_folder}_realigned${loc}.bam

#Index realignement
#samtools index ./mapping_PilotJan2015_realigned${loc}.bam


