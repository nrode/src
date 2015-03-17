#!/bin/bash

for file in *.fq.gz;do

filename=$(basename $file .fq.gz)
mv $file ${filename}.fastq.gz

done


