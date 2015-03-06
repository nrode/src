#!/bin/bash
#!/usr/local/bioinfo/R/bin/Rscript

#Aim: perform a fastqc analysis on a folder containing either the rawdata (argument1=1) or the demultiplexed data (arument1=2) in a given folder (argument 2)

if [ "$1" -ne "1" ] && [ "$1" -ne "2" ] && [ "$1" -ne "3" ] || [ -z $2 ]
then

echo "usage: \$1 =  1: analysis on raw_data (by  index)"
echo "or 2: analysis on demultiplexed data (by index/tag)"
echo "or 3: analysis on demultiplexed and trimmed data (by index/tag)echo "
echo "\$2 = name of the project folder"

exit 
fi

project_folder=$2
cd  ~/sead/projects/$project_folder
if (( $1==1 )); then
	echo "analysis on raw_data"
	qualitydir="qualityraw_${project_folder}"
	echo $qualitydir
	pathtodata=~/sead/raw_data/${project_folder}/
elif (( $1==2 )); then
	echo "analysis on demultiplexed data"
	qualitydir="qualitydemultiplex_${project_folder}"
	pathtodata=~/sead/projects/${project_folder}/demultiplex_${project_folder}/
else
	echo "analysis on trimmed data"
	qualitydir="qualityposttrim_${project_folder}"
	pathtodata="~/sead/projects/${project_folder}/trimmomatic_${project_folder}/"
fi

if [[ -d $qualitydir ]]; then
	echo "Former folder ($qualitydir) has been removed."
	rm -r $qualitydir;
fi

mkdir $qualitydir
cd $qualitydir


index=$(awk '{print $1}' ../tags/index.txt)

#index="Index_10"

for i in $index; do
		seqindex=$(awk -v row=$i '{if($1==row){ print $2 }}' ../tags/index.txt)
		#R2 adaptor sequence on nanuq
		echo -e "PE1rc\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" >> adapter.fasta
		#R1 adaptor sequence on nanuq
		echo -e "PE2rc\t"${seqindex} >> adapter.fasta
		#Adapter sequences from Sylvain
		#echo -e "PE1\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> adapter.fasta
		#echo -e "PE2\tCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> adapter.fasta
		#echo -e "${j}PE1rc\t"${seqtagrev}"AGATCGGAAGAGCACACGTCTGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTA" >> adapter.fasta
		#echo -e "${j}PE2rc\t"${seqtagrev}"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" >> adapter.fasta

		echo -e "Adaptater sequences for forward and reverse reads (tag: "${i}") have been created."
		
		if(($1==1)); then

			for k in ${pathtodata}*${i}*.fastq.gz;do 
				
				#Analysis with 13bp Illumina adapters: AGATCGGAAGAGC
				#/usr/local/bioinfo/FastQC/fastqc -o ./ -t 2 $k
				#Analysis with user-provided adapters
				/usr/local/bioinfo/FastQC/fastqc -o ./ -a adapter.fasta -t 2 $k
			done
			rm adapter.fasta
		fi

done

echo "fastqc done"

#Merge png quality files
../../../../src/quantgen/fastqc_cat-png.bash -I "*_fastqc.zip" -p per_base_quality -o plots_per_base_quality.pdf -c

cd ../../../src

echo $qualitydir
Rscript summary_fastqc.R $project_folder $qualitydir

if (( $1==3 )); then
cd ~/sead/projects/${project_folder}/$qualitydir
join -1 1 -2 1 ../qualitydemultiplex_${project_folder}/readcounts.txt ./readcounts.txt > summaryreadcountsposttrim.txt

fi
