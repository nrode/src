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

with V1.0 use -a AGATCGGAAGAG -A TGGAAAGATCGGAAGAG -o bank1/0Mtp1005_clean_R1.fq.gz -p bank1/0Mtp1005_clean_R2.fq.gz -e 0.1 -O 3 -m 35 -q 20 --max-n 0.2 ../demultiplexing_tim/bank1/Vigne-1_NoIndex_L005_m et4c_0Mtp1005_R1.fastq.gz ../demultiplexing_tim/bank1/Vigne-1_NoIndex_L005_met4c_0Mtp1005_R2.fastq .gz


