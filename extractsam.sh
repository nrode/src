#!/bin/bash

cd /home/rode/sead/projects/PilotJan2015/mapping_PilotJan2015
bamfile="mapping_PilotJan2015.bam"
#With regular expression and grep (fast: ~5min per sam)
samtools view $bamfile | grep "RPI[2345]" >> mapping_PilotJan2015RPI2345reggrep.sam
samtools view $bamfile | grep "RPI[6789]" >> mapping_PilotJan2015RPI6789reggrep.sam

#With egrep (fast: ~6min per sam)
#samtools view $bamfile | egrep "RPI2|RPI3|RPI4|RPI5" >> mapping_PilotJan2015RPI2345grep.sam
#samtools view $bamfile | egrep "RPI6|RPI7|RPI8|RPI9" >> mapping_PilotJan2015RPI6789grep.sam

#With awk pb since all the read groups are not in the same column
#samtools view $bamfile | awk '{if($13 ~ "RPI2" || $13 ~ "RPI3" || $13 ~ "RPI4" || $13 ~ "RPI5"){print $0}}' >> testmapping_PilotJan2015RPI2345.sam
#samtools view $bamfile | awk '{if($13 ~ "RPI6" || $13 ~ "RPI7" || $13 ~ "RPI8" || $13 ~ "RPI9"){print $0}}' >> mapping_PilotJan2015RPI6789.sam


#With awk match anywhere in on a given raw (slow: ~16 min per sam)
#samtools view $bamfile | awk '{if(/RPI2/ || /RPI3/ || /RPI4/ || /RPI5/){print $0}}' >> mapping_PilotJan2015RPI2345.sam
#samtools view $bamfile | awk '{if(/RPI6/ || /RPI7/ || /RPI8/ || /RPI9/){print $0}}' >> mapping_PilotJan2015RPI6789.sam

head -n 1000000 mapping_PilotJan2015RPI2345.sam | awk '{if($2==99 && $9<120){print $0}}'

