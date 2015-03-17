#!/bin/bash

#Aim: Extract all the information for a given locus from each bam file in a folder

project_folder=$1

pathtobam="../mapping_${project_folder}_indmerge"

loc="chr1"

cd ~/sead/projects/${project_folder}/mapping_${project_folder}/

dir=mapping_${project_folder}_$loc

if [[ -d $dir ]]; then
        echo "Former folder ($dir) has been removed."
        rm -r $dir;
fi

mkdir $dir
cd $dir



for bamfile in $pathtobam/*.bam;do

qsub -N split -cwd -b Y -q bioinfo.q ../../../../src/bamsplitloc.sh $pathtobam/$bamfile $loc

done

