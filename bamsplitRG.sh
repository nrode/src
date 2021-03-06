#!/usr/bin/env bash

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
  msg="\'${0##*/}' extracts a subset of readgroups from a bam file.\n"
  msg+="\n"
  msg+="Usage: ${0##*/} [OPTIONS] ...\n"
  msg+="\n"
  msg+="Options:\n"
  msg+="  -h, --help\tdisplay the help and exit\n"
  msg+="  -b, --bam\tname of the input bam file\n"
  msg+="  -i, --index\tname of the index/indices to extract in the bam\n"
  msg+="  -l, --locus\tname of the locus to extract in the bam (optional)\n"
  msg+="  -o, --out\toutput a list with the names of the bam file(s) created\n"
  msg+="\n"
  msg+="Examples:\n"
  msg+="  ${0##*/} -i RPI2 -t L1107 -o adapter.fasta -k -f\n"
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
          TEMP=`getopt -o hb:i:l:o -l help,bam:,index:,locus:,out, \
        -n "$0" -- "$@"`
  else # original getopt is available (no long options, whitespace, sorting)
          TEMP=`getopt hb:i:l:o "$@"`
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
      -b | --bam) bam=$2; shift 2;;
      -i | --index) index=$2; shift 2;;
      -l | --locus) locus=$2; shift 2;;
      -o | --out) out=true; shift;;
      --) shift; break;;
      *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
    esac
  done

  if [ -z "${index}" -a -z "${bam}" ]; then
    echo -e "ERROR: missing compulsory option -b or -i\n" 1>&2
    help
    exit 1
  fi

}

bam=""
index=""
locus=""
out=false

parseCmdLine "$@"

#Split a bam file according to the different read groups included in the header


bamname=$(basename $bam .bam)


#index="RPI2|RPI3|RPI5"
#index="RPI6|RPI7|RPI9"

#loc="scaffold0002"

#out=${4:-no}

if [[ "${out}" == true ]];then

	if [[ -e $out ]];then

		rm $out
		echo "Former list ("${out}") has been removed"
	fi
fi

RG=$(samtools view -H $bam | grep "^@RG" |cut -f2|cut -f2 -d':'| egrep $index)

echo $RG > RG.list
set -xv

ind=$(echo $index | tr '\|' '_')

filename=${bamname}_${ind}.bam

samtools view -h -b -R RG.list $bam $locus > $filename

	if [[ "${out}" == true ]];then

		echo $filename >> $out

	fi

rm RG.list
