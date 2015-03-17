#!/bin/bash
#!/usr/local/bioinfo/R/bin/Rscript
#Usage: Return one pdf file with the "per base sequence quality" and one pdf with the number of sequences for each original fastq files, the Phred/sequence length distributions and the percentage of adapter content.
#$1: name of the project
#$2: name of the directory with the fastqc.zip outputs

project_folder=$1
qualitydir=$2

#Merge png quality files
../../../../src/quantgen/fastqc_cat-png.bash -I "*_fastqc.zip" -p per_base_quality -o plots_per_base_quality.pdf -c

echo $qualitydir

. /etc/profile.d/modules.sh
module load compiler/gcc-4.8.2

Rscript ../../../src/summary_fastqc.R $project_folder $qualitydir


