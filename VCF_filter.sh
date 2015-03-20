#!/bin/bash

project_folder=$1

cd ~/sead/projects/${project_folder}

pathtovcf=" /home/rode/sead/projects/PilotJan2015/haplocall_PilotJan2015chr1"

dir="haplofilter_${project_folder}chr1"

if [[ -d $dir ]]; then
        echo "Former folder ($dir) has been removed."
        rm -r $dir;
fi



mkdir $dir

cd $dir

for file in $pathtovcf/*.vcf;do

filename=$(basename $file .vcf)

qsub -b Y -N VCFfilter -cwd -q bioinfo.q arcad_hts_Filter_VCF_on_individual_depth.pl -i $file -o ./${filename}_filtered.vcf -ind_prefix RPI -ind_min_depth 100 -ind_min_number 14

done

