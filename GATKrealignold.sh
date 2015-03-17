#!/bin/bash


#samtools view -f 4 mapping_PilotJan2015.bam chr1 | wc -l

#samtools idxstats mapping_PilotJan2015.bam | head

#samtools flagstat mapping_PilotJan2015.bam


#perl /home/sarah1/src/arcad-hts/arcad-hts/scripts/arcad_hts_6_bamRealignerRecalibrate.pl -h


/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator --help

#With encoding option (doesn't work)
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -o realign_PilotJan2015/real/RTC.intervals -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -fixMisencodedQuals -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta

#Look for positions to be realigned
#Without encoding option
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -o realign_PilotJan2015/real/RTC.intervals -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta


#Perform the realignment
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -targetIntervals realign_PilotJan2015/real/RTC.intervals -o realign_PilotJan2015/mapping_PilotJan2015_realigned.bam

#Index realignement
samtools index realign_PilotJan2015/mapping_PilotJan2015_realigned.bam


#cd SNPCallingbcftools

#SNP Calling with a pipe
#samtools mpileup -ugf ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta ../realign_PilotJan2015/mapping_PilotJan2015_realigned.bam | bcftools call -vmO z -o PilotJan2015.vcf.gz

#a/SNPCalling in two steps
#samtools mpileup -go PilotJan2015.bcf -f ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta ../realign_PilotJan2015/mapping_PilotJan2015_realigned.bam
#bcftools call -vmO z -o PilotJan2015.vcf.gz PilotJan2015.bcf

#File preparation
tabix -p vcf PilotJan2015.vcf.gz

#Making graph to filter the variants
bcftools stats -F ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -s - PilotJan2015.vcf.gz Pilotjan2015.vcf.gz.stats
mkdir plots
plot-vcfstats -p plots/Pilotjan2015.vcf.gz.stats

#Filter the data
bcftools filter -O z -o PilotJan2015_filtered..vcf.gz -s LOWQUAL -i'%QUAL>10' PilotJan2015.vcf.gz

#b/SNP calling using haplotype caller
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ../realign_PilotJan2015/mapping_PilotJan2015_realigned.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o raw_variants.vcf






cd realign_PilotJan2015
mkdir recal
chmod 775

cd ..
#Prepare table for base recalibration
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -knownSites dbSNP -I realign_PilotJan2015/mapping_PilotJan2015_realigned.bam -o realign_PilotJan2015/recal/mapping_PilotJan2015_recal.table

#Prepare table after first recalibration
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I realign_PilotJan2015/mapping_PilotJan2015_realigned.bam --BQSR realign_PilotJan2015/recal/mapping_PilotJan2015_recal.table -o realign_PilotJan2015/recal/mapping_PilotJan2015_post_recal.table 

#Check base recalibration using install.packages("gsalib") in  R studio
#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T AnalyzeCovariates -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I realign_PilotJan2015/mapping_PilotJan2015_realigned.bam -before realign_PilotJan2015/recal/mapping_PilotJan2015_recal.table -after realign_PilotJan2015/recal/mapping_PilotJan2015_post_recal.table -plots realign_PilotJan2015/recal/recalibration_plots.pdf



#Perform base recalibration
/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T PrintReads -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I realign_PilotJan2015/mapping_PilotJan2015_realigned.bam --BSQR realign_PilotJan2015/recal/mapping_PilotJan2015_recal.table -o realign_PilotJan2015/mapping_PilotJan2015_recal.bam




#Realign using the Arcard script (doesn't work)
qsub -N realign -b Y -cwd -q bioinfo.q perl /home/sarah1/src/arcad-hts/arcad-hts/scripts/arcad_hts_6_bamRealignerRecalibrate.pl -I ./mapping_PilotJan2015/mapping_PilotJan2015.bam -R ../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -o realign_PilotJan2015/



#all SNPs and short indels for multiple diploid individuals:

#samtools mpileup -P ILLUMINA -ugf ref.fa *.bam | bcftools view -bcvg - > var.raw.bcf
#bcftools view var.raw.bcf | vcfutils.pl varFilter -D 2000 > var.flt.vcf

#Individuals are identified from the SM tags in the @RG header lines. Individuals can be pooled in one alignment file; one individual can also be separated into multiple files. The -P option specifies that indel candidates should be collected only from read groups with the @RG-PL tag set to ILLUMINA. Collecting indel candidates from reads sequenced by an indel-prone technology may affect the performance of indel calling. 


