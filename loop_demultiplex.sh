#!/bin/bash
# usage: ~/rode/sead/src/loop_demultiplex.sh project_folder  
# arguments: $1 = name of the project's folder
# tags descriptions (.txt files) should be in the "tags" folder in "project_folder"

project_folder=$1
cd ~/sead/projects/$project_folder
if [[ -d demultiplex_${project_folder} ]]; then
	echo "Former folder demultiplex_${project_folder} has been removed."
	rm -r demultiplex_${project_folder};
fi
mkdir demultiplex_${project_folder}
cd demultiplex_${project_folder}
s=0
for f in ../tags/forward_RPI*.txt; do
	((s=s+1))
	index=$(basename $f .txt)
	index=${index#*_}
	if ((s<5)); then
		/home/rode/src/quantgen/demultiplex.py --idir ~/sead/raw_data/${project_folder}/ --ifq1 ${index}_L001_R1_001.fastq.gz --ifq2 ${index}_L001_R2_001.fastq.gz --it $f -v 1 --ofqp $index --met 4d --dist 20 --re ApeKI >> stdoutdemultiplex.txt
	else
		/home/rode/src/quantgen/demultiplex.py --idir ~/sead/raw_data/${project_folder}/ --ifq1 ${index}_L001_R1_001.fastq.gz --ifq2 ${index}_L001_R2_001.fastq.gz --it $f -v 1 --ofqp $index --met 4d --dist 20 --re EcoT22I >> stdoutdemultiplex.txt
	fi
done
