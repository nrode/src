#!/bin/bash

project_folder=$1

cd ~/sead/projects/${project_folder}

mapdir="mapping_${project_folder}"

if [[ -d $mapdir ]]; then
        echo "Former folder ($mapdir) has been removed."
        rm -r $mapdir;
fi

mkdir $mapdir
cd trimmomatic_${project_folder}

for i in *R1.paired.fastq.gz; do
	rev=${i/R1/R2}
	here=$(pwd)
	name=${i%_*}
	#pairing=${i#*.}

	echo -e "$here/$i\t$name\tpaired\t$here/$rev" >> conf_mapping.txt

	unpairedR1=${i/paired.fastq.gz/unpaired.fastq.gz}
	echo -e "$here/$unpairedR1\t$name\tsingle" >> conf_mapping.txt
	
	unpairedR2=${rev/paired.fastq.gz/unpaired.fastq.gz}
	echo -e "$here/$unpairedR2\t$name\tsingle" >> conf_mapping.txt

done

mv ./conf_mapping.txt ../$mapdir
cd ../$mapdir 


perl /home/sarah1/src/arcad-hts/arcad-hts/scripts/arcad_hts_4_Mapping_Arcad.pl -i conf_mapping.txt -o mapping_${project_folder}.bam -r /home/rode/sead/raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -q bioinfo.q -normdup -mapper bwa_mem

perl /home/sarah1/src/arcad-hts/arcad-hts/scripts/arcad_hts_5_cleanMapping.pl --bam /home/rode/sead/projects/${project_folder}/mapping_${project_folder}/mapping_${project_folder}.bam

