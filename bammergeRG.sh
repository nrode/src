#!/bin/bash

project_folder=$1

cd ~/sead/projects/${project_folder}

dir="mapping_${project_folder}_indmerge"

if [[ -d $dir ]]; then
        echo "Former folder ($dir) has been removed."
        rm -r $dir;
fi

cd mapping_${project_folder}/mapping_${project_folder}_tmp

mkdir ../$dir

for file in *.paired.bam;do

name=$(basename $file .paired.bam)

samtools merge ../$dir/$name.bam $file $name.single.bam $name.single_1.bam

done


