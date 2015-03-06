#!/bin/bash

samfile=$1
motif=$2
#flag=${3:-""}


#samtools view $bamfile | awk '{if($13 ~ "RPI2" || $13 ~ "RPI3" || $13 ~ "RPI4" || $13 ~ "RPI5"){print $0}}' >> mapping_PilotJan2015RPI2345.sam
#samtools view $bamfile | awk '{if($13 ~ "RPI6" || $13 ~ "RPI7" || $13 ~ "RPI8" || $13 ~ "RPI9"){print $0}}' >> mapping_PilotJan2015RPI6789.sam

#tot=$(samtools view $flag $bamfile | grep -c "." )
#num=$(samtools view $flag $bamfile | grep -c $motif )

tot=$(grep -c "." $samfile)
num=$(grep -c $motif $samfile)
perc=$(echo "$num/$tot" | bc -l) 

echo $samfile $motif $num $tot $perc
