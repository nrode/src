# src
Scripts for Genotyping By Sequencing (GBS) Analyses


chimeracount.sh: count the number of sequences with a given motif in a fastq file

createadapter.sh: recreate the sequences of custom-made adapters in a file to be passed on to other programs such as fasqc (default text format) or trimmomatic (option -f: fasta format)

loop.fastqcadaptnew.sh: fasqc analysis of of different fastqc files/compute summary graphs using the script summary_fastqc.R

summary_fastqc.R: script combining different functions from other github scripts (https://github.com/timflutre/quantgen/blob/master/fastqc_cat-png.bash and https://github.com/timflutre/quantgen/blob/master/utils_fastqc.R)

