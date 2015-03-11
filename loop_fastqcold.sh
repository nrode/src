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
	datatoanalyse=~/sead/raw_data/${project_folder}/RPI2_*.fastq.gz
elif (( $1==2 )); then
	qualitydir="qualitydemultiplex_${project_folder}"
	datatoanalyse=~/sead/projects/${project_folder}/demultiplex_${project_folder}/*.fastq.gz  
else
	qualitydir="qualityposttrim_${project_folder}"
	datatoanalyse=~/sead/projects/${project_folder}/trimmomatic_${project_folder}/*.fastq.gz
fi
if [[ -d $qualitydir ]]; then
	echo "Former folder ($qualitydir) has been removed."
	rm -r $qualitydir;
fi
echo $qualitydir
mkdir $qualitydir
cd $qualitydir
/usr/local/bioinfo/FastQC/fastqc -o ./ -t 2 $datatoanalyse
echo "fastqc done"

#Merge png quality files
../../../../src/quantgen/fastqc_cat-png.bash -I "*_fastqc.zip" -p per_base_quality -o plots_per_base_quality.pdf -c


cd ../../../src
echo $qualitydir
pwd

Rscript summary_fastqc.R $project_folder $qualitydir

if (( $1==3 )); then
cd ~/sead/projects/${project_folder}/$qualitydir
join -1 1 -2 1 ../qualitydemultiplex_${project_folder}/readcounts.txt ./readcounts.txt > summaryreadcountsposttrim.txt

fi
