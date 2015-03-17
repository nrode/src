

/usr/local/jre/bin/java -jar /usr/local/bioinfo/Trimmomatic/trimmomatic-0.33.jar PE -trimlog trim.txt test.fastq test2.fastq testpairedR1.fq.gz testpairedR2.fq.gz testunpairedR1.fq.gz testunpairedR2.fq.gz ILLUMINACLIP:adapter.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
