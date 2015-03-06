#!/bin/bash

project_folder=$1

cd /home/rode/sead/projects/${project_folder}

cleanmapdir="cleanmapping_${project_folder}"

if [[ -d $cleanmapdir ]]; then
        echo "Former folder ($cleanmapdir) has been removed."
        rm -r $cleanmapdir;
fi

echo ${cleanmapdir}
mkdir ${cleanmapdir}
cd ~/sead/projects/${project_folder}/${cleanmapdir}
#mypath="~/sead/projects/${project_folder}/${cleanmapdir}"

perl /home/sarah1/src/arcad-hts/arcad-hts/sp1_snp_detection_on_rnaseq/5_cleanMapping.pl --bam /home/rode/sead/projects/PilotJan2015/mapping_${project_folder}/mapping_${project_folder}.bam
