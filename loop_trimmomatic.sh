#!/bin/bash

project_folder=$1

cd ~/sead/projects/${project_folder}

trimdir="trimmomatic_${project_folder}"

if [[ -d $trimdir ]]; then
        echo "Former folder ($trimdir) has been removed."
        rm -r $trimdir;
fi

mkdir $trimdir
cd $trimdir

echo -e "indextag\tinitialreadnumber\ttrimfromstart\ttrimfromend\ttooshort\treadnumberaftertrim\tpairedreadnumber\tunpairedR1readnumber\tunpairedR2readnumber" > summarytrim.txt

index=$(awk '{print $1}' ../tags/index.txt)

#index="RPI2"

for i in $index; do
((s=$(expr substr $i 4 4)-1))

mypath="../tags/${i}_*.txt"
tag=$(awk '{print $1}' $mypath)
#tag="L1107"

	seqindex=$(awk -v row=$i '{if($1==row){ print $2 }}' ../tags/index.txt)
	seqindexrev=$(awk -v row=$i '{if($1==row){ print $3 }}' ../tags/index.txt)
	echo $seqindex
	echo $seqindexrev

	for j in $tag;do
		seqtag=$(awk -v row=$j '{if($1==row){ print $2 }}' ../tags/${i}_*.txt)
		echo $seqtag
		seqtagrev=$(awk -v row=$j '{if($1==row){ print $3 }}' ../tags/${i}_*.txt)
		indextag=${i}_S${s}_${j}
		echo -e ">Prefix_PE/1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} > adapter.fasta
		echo -e ">Prefix_PE/2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> adapter.fasta
		echo -e ">PE1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> adapter.fasta
		echo -e ">PE2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> adapter.fasta
		echo -e ">PE1rc\n"${seqtagrev}"AGATCGGAAGAGCACACGTCTGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTA" >> adapter.fasta
		echo -e ">PE2rc\n"${seqtagrev}"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" >> adapter.fasta

		echo -e "Adaptater sequences for forward and reverse reads (tag: "${j}") have been created."
	
		echo "../demultiplex_${project_folder}/${indextag}_R1.fastq.gz "
		echo "../demultiplex_${project_folder}/${indextag}_R2.fastq.gz "
		
		((minlength=36))		
		/usr/local/jre/bin/java -jar /usr/local/bioinfo/Trimmomatic/trimmomatic-0.33.jar PE -trimlog ${indextag}_trimlog.txt ~/sead/projects/${project_folder}/demultiplex_${project_folder}/${indextag}_R1.fastq.gz ~/sead/projects/${project_folder}/demultiplex_${project_folder}/${indextag}_R2.fastq.gz ${indextag}_R1.paired.fastq.gz ${indextag}_R1.unpaired.fastq.gz ${indextag}_R2.paired.fastq.gz ${indextag}_R2.unpaired.fastq.gz ILLUMINACLIP:adapter.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minlength
		
		readnumber=$(wc -l ${indextag}_trimlog.txt | cut -d " " -f1)
		trimfromstart=$(awk 'BEGIN {s=0}{if($4>0){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		trimfromend=$(awk 'BEGIN {s=0}{if($6>0){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		tooshort=$(awk -v l=$minlength 'BEGIN {s=0}{if($3<l){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		((readnumberaftertrim=$readnumber-$tooshort))
		pairedreadnumber=$(zcat ${indextag}_R1.paired.fastq.gz |awk ' END {print NR/4}')
		unpairedR1readnumber=$(zcat ${indextag}_R1.unpaired.fastq.gz |awk ' END {print NR/4}')
		unpairedR2readnumber=$(zcat ${indextag}_R2.unpaired.fastq.gz |awk ' END {print NR/4}')


		echo -e ${indextag}"\t"${readnumber}"\t"${trimfromstart}"\t"${trimfromend}"\t"${tooshort}"\t"${readnumberaftertrim}"\t"${pairedreadnumber}"\t"${unpairedR1readnumber}"\t"${unpairedR2readnumber} >> summarytrim.txt	

		done

done



