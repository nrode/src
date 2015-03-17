#!/bin/bash
#!/usr/local/bioinfo/R/bin/Rscript
#Aim: compare the number of base pairs in paired-end fastq files in a folder containing either the rawdata (argument1=1) or the demultiplexed data (arument1=2) in a given folder (argument 2)

if [ "$1" -ne "1" ] && [ "$1" -ne "2" ] && [ "$1" -ne "3" ] || [ -z $2 ]
then

echo "usage: \$1 =  1: analysis on raw_data"
echo "or 2: analysis on demultiplexed data"
echo "or 3: analysis on demultiplexed and trimmed data"
echo "\$2 = name of the project folder"

exit 
fi

project_folder=$2
cd  ~/sead/projects/${project_folder}
if (( $1==1 )); then
	echo "analysis on raw_data"
	dir="seqlenraw_${project_folder}"
	echo $dir
	pathtodata=../../../raw_data/${project_folder}
elif (( $1==2 )); then
	echo "analysis on demultiplexed data"
	dir="seqlendemultiplex_${project_folder}"
	pathtodata=../demultiplex_${project_folder}
else
	echo "analysis on trimmed data"
	dir="seqlenposttrim_${project_folder}"
	pathtodata=../trimmomatic_${project_folder}
fi

if [[ -d $dir ]]; then
	echo "Former folder ($dir) has been removed."
	rm -r $dir
fi

mkdir $dir
cd $dir
set -xv



#name="*_R1*fastq.gz"

#for file in ${name};do
for file in RPI*_L001_R1_001.fastq.gz;do
#for file in RPI6_S5_L1107_R1.fastq.gz;do

	qsub -V -b Y -cwd -q bioinfo.q "../../../src/lengthseqfastq.sh $pathtodata/$file"
	#qsub sleep 100 -N $file -b Y -cwd -q bioinfo.q "../../../src/lengthseqfastq.sh $pathtodata/$file"


done

#Histogram with R
#. /etc/profile.d/modules.sh
#module load compiler/gcc-4.8.2
#qsub -N "HoldJob" -hold_jid "RPI*" -b y -q bioinfo.q -cwd Rscript ../../../src/summary_lengthfastq.R $project_folder $dir



#qsub sleep 30 -N "HoldJob" -hold_jid "RPI*" -b y -q bioinfo.q -cwd "Rscript ../../../src/summary_lengthfastq.R "

