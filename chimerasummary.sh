#!/bin/bash

cd /home/rode/sead/projects/PilotJan2015/mapping_PilotJan2015

samfile=mapping_PilotJan2015RPI2345grep.sam

/home/rode/sead/src/chimeracount.sh $samfile "GC[AT]GC" >> chimeracount.out
/home/rode/sead/src/chimeracount.sh $samfile "GCAGC"  >> chimeracount.out
/home/rode/sead/src/chimeracount.sh $samfile "GCTGC"  >> chimeracount.out


samfile=mapping_PilotJan2015RPI6789grep.sam
/home/rode/sead/src/chimeracount.sh $samfile "ATGCAT"  >> chimeracount.out

samfile=mapping_PilotJan2015RPI2345grepproper.sam

/home/rode/sead/src/chimeracount.sh $samfile "GC[AT]GC" >> chimeracount.out
/home/rode/sead/src/chimeracount.sh $samfile "GCAGC"  >> chimeracount.out
/home/rode/sead/src/chimeracount.sh $samfile "GCTGC"  >> chimeracount.out


samfile=mapping_PilotJan2015RPI6789grepproper.sam
/home/rode/sead/src/chimeracount.sh $samfile "ATGCAT"  >> chimeracount.out


