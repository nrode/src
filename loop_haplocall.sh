#!/usr/bin/env bash

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
  msg="\'${0##*/}' make a loop over all the bam files in an input folder and call haplotypes using Haplotype Caller.\n"
  msg+="\n"
  msg+="Usage: ${0##*/} [OPTIONS] ...\n"
  msg+="\n"
  msg+="Options:\n"
  msg+="  -h, --help\tdisplay the help and exit\n"
  msg+="  -f, --folder\tname of the folder with the different bam files\n"
  msg+="  -s, --suffix\tname of the suffix to be appended at the end of each vcf file created (default = raw)\n"
  msg+="\n"
  msg+="Examples:\n"
  msg+="  ${0##*/} -f PilotJan2015 -o raw.variant.vcf \n"
  msg+="\n"
  msg+="Report bugs to <nicolas.o.rode@gmail.com>."
  echo -e "$msg"
}

# Parse the command-line arguments.
# http://stackoverflow.com/a/4300224/597069
# http://stackoverflow.com/questions/402377/


function parseCmdLine () {
  getopt -T > /dev/null # portability check (say, Linux or Mac OS?)
  if [ $? -eq 4 ]; then # GNU enhanced getopt is available
          TEMP=`getopt -o hf:s -l help,index:,keep,fasta, \
        -n "$0" -- "$@"`
  else # original getopt is available (no long options, whitespace, sorting)
          TEMP=`getopt hf:s "$@"`
  fi
  if [ $? -ne 0 ]; then
          echo "ERROR: "$(which getopt)" failed" 1>&2
          getopt -T > /dev/null
          if [ $? -ne 4 ]; then
            echo "did you use long options? they are not handled \
on your system, use -h for help"
          fi
          exit 2
  fi
  eval set -- "$TEMP"
  while [ $# -gt 0 ]; do
    case "$1" in
      -h | --help) help; exit 0; shift;;
      -f | --folder) folder=$2; shift 2;;
      -s | --suffix) suffix=$2; shift 2;;
      --) shift; break;;
      *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
    esac
  done

  if [ -z "${folder}" ]; then
    echo -e "ERROR: missing compulsory option -f \n" 1>&2
    help
    exit 1
  fi

}

folder=""
suffix="raw"

parseCmdLine "$@"

project_folder=$folder

cd ~/sead/projects/${project_folder}

dir="haplocall_${project_folder}"

if [[ -d $dir ]]; then
        echo "Former folder ($dir) has been removed."
        rm -r $dir;
fi

mkdir $dir

cd $dir

pathtobam=../bamsplit_${project_folder}
inbam=mapping_${project_folder}_realigned.bam

for input in $pathtobam/*.bam;do

echo $input
inputname=$(basename $input .bam)
#SNP calling using haplotype caller

#/usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I $pathtobam/${inbam} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -min_mapping_quality_score 30 -o raw_variants.vcf

samtools index $input

qsub -b Y -cwd -q bioinfo.q -N haplocall /usr/local/jre/bin/java -jar /usr/local/bioinfo/GATK/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../../raw_data/RefGenome/JCVI.Medtr.v4.20130313.fasta -I ${input} --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ${inputname}${suffix}.vcf


done
