#!/bin/bash

project_folder=$1

cd ~/sead/projects/${project_folder}/cutadapt_${project_folder}
ls
index=$(awk '{print $1}' index.txt)
#index="RPI2"

for i in $index; do

if [[ -d $i ]]; then
	echo "Former folder ($i) has been removed."
	rm -r $i;
fi

mkdir $i
cd $i

seqindex=$(awk -v row=$i '{if($1==row){ print $2 }}' ../index.txt)

echo -e ">adapter\nAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" > adapter.fasta;
echo -e "Adaptater sequence (tag: "${seqindex}") has been created."

path=../../demultiplex_${project_folder}/

	for j in ${path}/${i}*.fastq.gz; do
	filename=$(basename ${j} .gz)	
	gunzip -c ${j} > $filename	
	done
cd ..
done


perl /home/sarah1/src/arcad-hts/arcad-hts/scripts/arcad_hts_1_cutadapt_in_chain.pl -i . -sub -sb adapter.fasta -o posttrim


