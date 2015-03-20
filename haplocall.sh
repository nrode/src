#!/bin/bash
project_folder=$1


cd ~/sead/projects/${project_folder}

dir="haplocall_${project_folder}_chr1"

if [[ -d $dir ]]; then
        echo "Former folder ($dir) has been removed."
        rm -r $dir;
fi

mkdir $dir

cd $dir

pathtobam=../realign_${project_folder}chr1
input=mapping_${project_folder}_realigned.bam
#SNP calling using haplotype caller

#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I $pathtobam/${input} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -min_mapping_quality_score 30 -o raw_variants.vcf


/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I $pathtobam/${input} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -min_mapping_quality_score 30 -o raw_variants.vcf



