#!/bin/bash
#!/usr/local/bioinfo/R/bin/Rscript
set -xv
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
cd  ~/sead/projects/${project_folder}
if (( $1==1 )); then
	echo "analysis on raw_data"
	qualitydir="qualityraw_${project_folder}"
	echo $qualitydir
	pathtodata=../../../raw_data/${project_folder}
elif (( $1==2 )); then
	echo "analysis on demultiplexed data"
	qualitydir="qualitydemultiplex_${project_folder}"
	pathtodata=./demultiplex_${project_folder}/
else
	echo "analysis on trimmed data"
	qualitydir="qualityposttrim_${project_folder}"
	pathtodata=./trimmomatic_${project_folder}/
fi

if [[ -d $qualitydir ]]; then
	echo "Former folder ($qualitydir) has been removed."
	rm -r $qualitydir;
fi

mkdir $qualitydir
cd $qualitydir


#index=$(awk '{print $1}' ../tags/index.txt)

index="RPI2"

#Loop over each index
for indexname in $index; do
	
	s=$(($(expr substr $indexname 4 4)-1))

	pathtotags="../tags/${indexname}_*.txt"
	tag=$(awk '{print $1}' $pathtotags)
	#tag="L1107"

	#Add Illumina Universal Adapter at the begining of the adapter file
	echo -e "Universal_Illumina_Adapter\tAGATCGGAAGAG" >> adapter.txt

	
	for tagname in $tag;do

		indextag=${indexname}_S${s}_${tagname}

		if (($1==1));then
	
			#Keep previous file and append sequences for this index/tag combination
			../../../src/createadapter.sh -i $indexname -t $tagname -o adapter.txt -k
	
		else
			#Erase previous file and create a new file for this index/tag combination
			../../../src/createadapter.sh -i $indexname -t $tagname -o adaper.txt

			for file in ${pathtodata}/${indextag}_R*.fastq.gz;do 
				
				#Analysis with in-house adapters
				/usr/local/bioinfo/FastQC/fastqc -o ./ -a adapter.txt -t 2 $file

			done
	
		fi 
	done

	if (($1==1));then
		
		#Analysis with in-house adapters
		/usr/local/bioinfo/FastQC/fastqc -o ./ -a ./adapter.txt -t 2 ${pathtodata}/${indexname}*.fastq.gz
	fi

done

echo "fastqc done"

#Merge png quality files
../../../../src/quantgen/fastqc_cat-png.bash -I "*_fastqc.zip" -p per_base_quality -o plots_per_base_quality.pdf -c

echo $qualitydir
Rscript ../../../src/summary_fastqc.R $project_folder $qualitydir

if (( $1==3 )); then
cd ~/sead/projects/${project_folder}/$qualitydir
join -1 1 -2 1 ../qualitydemultiplex_${project_folder}/readcounts.txt ./readcounts.txt > summaryreadcountsposttrim.txt

fi
