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

#index="RPI6"

echo $index

for indexname in $index; do

	s=$(($(expr substr $indexname 4 4)-1))

	pathtotags="../tags/${indexname}_*.txt"
	tag=$(awk '{print $1}' $pathtotags)
	#tag="L1107"


	for tagname in $tag;do
		
		indextag=${indexname}_S${s}_${tagname}
		
		../../../src/createadapter.sh -i $indexname -t $tagname -o adapter.fasta -f
	
		((minlength=36))		
		/usr/local/jre/bin/java -jar /usr/local/bioinfo/Trimmomatic/trimmomatic-0.33.jar PE -trimlog ${indextag}_trimlog.txt ~/sead/projects/${project_folder}/demultiplex_${project_folder}/${indextag}_R1.fq.gz ~/sead/projects/${project_folder}/demultiplex_${project_folder}/${indextag}_R2.fq.gz ${indextag}_R1.paired.fq.gz ${indextag}_R1.unpaired.fq.gz ${indextag}_R2.paired.fq.gz ${indextag}_R2.unpaired.fq.gz ILLUMINACLIP:adapter.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:$minlength
		
		readnumber=$(wc -l ${indextag}_trimlog.txt | cut -d " " -f1)
		trimfromstart=$(awk 'BEGIN {s=0}{if($4>0){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		trimfromend=$(awk 'BEGIN {s=0}{if($6>0){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		tooshort=$(awk -v l=$minlength 'BEGIN {s=0}{if($3<l){s+=1}} END{print s}' ${indextag}_trimlog.txt)
		((readnumberaftertrim=$readnumber-$tooshort))
		pairedreadnumber=$(zcat ${indextag}_R1.paired.fq.gz |awk ' END {print NR/4}')
		unpairedR1readnumber=$(zcat ${indextag}_R1.unpaired.fq.gz |awk ' END {print NR/4}')
		unpairedR2readnumber=$(zcat ${indextag}_R2.unpaired.fq.gz |awk ' END {print NR/4}')


		echo -e ${indextag}"\t"${readnumber}"\t"${trimfromstart}"\t"${trimfromend}"\t"${tooshort}"\t"${readnumberaftertrim}"\t"${pairedreadnumber}"\t"${unpairedR1readnumber}"\t"${unpairedR2readnumber} >> summarytrim.txt	

		done

done



