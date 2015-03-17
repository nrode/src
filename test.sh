#!/bin/bash
#!/usr/local/bioinfo/R/bin/Rscript

. /etc/profile.d/modules.sh
module load compiler/gcc-4.8.2
Rscript ../../../src/summary_lengthfastq.R  PilotJan2015 seqlenraw_PilotJan2015
