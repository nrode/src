#!/bin/bash

project_folder=$1

pathtodata=~/sead/raw_data/${project_folder}/

cd ~/work/projects/${project_folder}

trimdir="trimmomatic_${project_folder}"

if [[ -d $trimdir ]]; then
        echo "Former folder ($trimdir) has been removed."
#        rm -r $trimdir;
fi

#mkdir $trimdir
cd $trimdir

echo -e "indextag\tinitialreadnumber\ttrimfromstart\ttrimfromend\ttooshort\treadnumberaftertrim\tpairedreadnumber\tunpairedR1readnumber\tunpairedR2readnumber" > summarytrim.txt

index=$(awk '{print $1}' ../tags/index.txt)

#index="Index_1"


for i in $index; do
		
		seqindex=$(awk -v row=$i '{if($1==row){ print $2 }}' ../tags/index.txt)
		seqindexrev=$(awk -v row=$i '{if($1==row){ print $3 }}' ../tags/index.txt)
		#R2 adaptor sequence on nanuq
		echo -e ">PrefixPE/1\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" > adapter.fasta
		#R1 adaptor sequence on nanuq
		echo -e ">PrefixPE/2\n"${seqindex} >> adapter.fasta
		#Other sequences to test in single mode
		echo -e ">PE1\nAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" >> adapter.fasta
		echo -e ">PE1rc\nACACTCTTTCCCTACACGACGCTCTTCCGATCT" >> adapter.fasta
		echo -e ">PE2\n"${seqindex} >> adapter.fasta
		echo -e ">PE2rc\n"${seqindexrev} >> adapter.fasta

		echo -e "Adaptater sequences for forward and reverse reads (index: "${i}") have been created."
		

			for k in ${pathtodata}*${i}.*_R1.fastq.gz;do 
				echo $k
				fullnameR1=$(basename $k .fastq.gz)
				fullnameR2=${fullnameR1/R1/R2}
				name=$(echo "Index_"$(echo $fullnameR1| cut -d'_' -f 2))
				nameR1=$(echo $name".R1")
				nameR2=${nameR1/R1/R2}

				#Analysis with user-provided adapters
				((minlength=36))		
			#	/usr/local/jre/bin/java -jar /usr/local/bioinfo/Trimmomatic/trimmomatic-0.33.jar PE -trimlog ${name}_trimlog.txt ~/sead/raw_data/${project_folder}/${fullnameR1} ~/sead/raw_data/${project_folder}/${fullnameR2} ${name}_R1.paired.fastq.gz ${name}_R1.unpaired.fastq.gz ${name}_R2.paired.fastq.gz ${name}_R2.unpaired.fastq.gz ILLUMINACLIP:adapter.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minlength

		readnumber=$(wc -l ${name}_trimlog.txt | cut -d " " -f1)
		trimfromstart=$(awk 'BEGIN {s=0}{if($4>0){s+=1}} END{print s}' ${name}_trimlog.txt)
		trimfromend=$(awk 'BEGIN {s=0}{if($6>0){s+=1}} END{print s}' ${name}_trimlog.txt)
		tooshort=$(awk -v l=$minlength 'BEGIN {s=0}{if($3<l){s+=1}} END{print s}' ${name}_trimlog.txt)
		((readnumberaftertrim=$readnumber-$tooshort))
		pairedreadnumber=$(zcat ${name}_R1.paired.fastq.gz |awk ' END {print NR/4}')
		unpairedR1readnumber=$(zcat ${name}_R1.unpaired.fastq.gz |awk ' END {print NR/4}')
		unpairedR2readnumber=$(zcat ${name}_R2.unpaired.fastq.gz |awk ' END {print NR/4}')

		echo -e ${name}"\t"${readnumber}"\t"${trimfromstart}"\t"${trimfromend}"\t"${tooshort}"\t"${readnumberaftertrim}"\t"${pairedreadnumber}"\t"${unpairedR1readnumber}"\t"${unpairedR2readnumber} >> summarytrim.txt	
				done

		rm adapter.fasta

done
