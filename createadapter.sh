#!/usr/bin/env bash

# Display the help on stdout.
# The format complies with help2man (http://www.gnu.org/s/help2man)
function help () {
  msg="\'${0##*/}' creates a sequence file with adapter sequences for an index/tag combination.\n"
  msg+="\n"
  msg+="Usage: ${0##*/} [OPTIONS] ...\n"
  msg+="\n"
  msg+="Options:\n"
  msg+="  -h, --help\tdisplay the help and exit\n"
  msg+="  -i, --index\tindex name of the sequence\n"
  msg+="  -t, --tag\ttag name of the sequence\n"
  msg+="  -o, --out\tname of the file created (default=adapter.fasta/adapter.txt)\n"
  msg+="  -k, --keep\tkeep a previous output file and append sequences to it\n"
  msg+="  -f, --fasta\tproduces a fasta instead of a text file\n"
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
          TEMP=`getopt -o hi:t:o:kf -l help,index:,tag:,out:,keep,fasta, \
        -n "$0" -- "$@"`
  else # original getopt is available (no long options, whitespace, sorting)
          TEMP=`getopt hi:t:o:kf "$@"`
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
      -i | --index) indexname=$2; shift 2;;
      -t | --tag) tagname=$2; shift 2;;
      -o | --out) outfile=$2; shift 2;;
      -k | --keep) keep=true; shift;;
      -f | --fasta) fasta=true; shift;;
      --) shift; break;;
      *) echo "ERROR: options parsing failed, use -h for help" 1>&2; exit 1;;
    esac
  done

  if [ -z "${indexname}" -a -z "${tagname}" ]; then
    echo -e "ERROR: missing compulsory option -i or -t\n" 1>&2
    help
    exit 1
  fi

  if [ -z "${outfile}" ]; then
    echo -e "ERROR: missing compulsory option -o\n" 1>&2
    help
    exit 1
  fi

}

indexname=""
tagname=""
outfile=""
keep=false
fasta=false

parseCmdLine "$@"

pathtoindexfile=../tags/index.txt
pathtotagfile=../tags/${indexname}_*.txt

	#Extract index sequence
	seqindex=$(awk -v row=$indexname '{if($1==row){ print $2 }}' $pathtoindexfile)
	seqindexrev=$(awk -v row=$indexname '{if($1==row){ print $3 }}' $pathtoindexfile)
	#Extract tag sequence
	seqtag=$(awk -v row=$tagname '{if($1==row){ print $2 }}' $pathtotagfile)
	seqtagrev=$(awk -v row=$tagname '{if($1==row){ print $3 }}' $pathtotagfile)
#	indextag=${i}_S${s}_${j}

if [[ -e $outfile ]];then
	
	if [[ "${keep}" == true ]];then

		echo "New sequences will be appended to the former file \""${outfile}"\"."

	else

		echo "Former file \""${outfile}"\" has been removed."
		rm $outfile
	fi

else
	echo "File \""${outfile}"\" has been created."
fi

if [[ "${fasta}" == true ]];then

	echo -e ">Prefix_PE/1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> $outfile
	echo -e ">Prefix_PE/2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> $outfile
	#echo -e ">PE1\nAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> $outfile
	#echo -e ">PE2\nCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> $outfile
	echo -e ">PE1rc\n"${seqtagrev}"AGATCGGAAGAGCACACGTCTGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTA" >> $outfile
	echo -e ">PE2rc\n"${seqtagrev}"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" >> $outfile

else
	if [[ ! -e $outfile ]];then

	echo -e "Universal_Illumina_Adapter\tAGATCGGAAGAG" >> $outfile

	fi

	#echo -e "PE1\tAATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"${seqtag} >> $outfile
	#echo -e "PE2\tCAAGCAGAAGACGGCATACGAGAT"${seqindexrev}"GTGACTGGAGTTCAGACGTGTCGTCTTCCGATCT"${seqtag} >> $outfile
	echo -e "TagRC6bp\t"${seqtagrev}"AGATCG" >> $outfile
	#echo -e "PE1_2rc${seqtagrev}+4bp\t"${seqtagrev}"AGAT" >> $outfile
	#echo -e "PE1_2rc${seqtagrev}+8bp\t"${seqtagrev}"AGATCGGA" >> $outfile
	#echo -e "PE1rc\t"${seqtagrev}"AGATCGGAAGAGCACACGTCTGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTA" >> $outfile
	#echo -e "PE2rc\t"${seqtagrev}"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"${seqindex}"ATCTCGTATGCCGTCTTCTCGTTGA" >> $outfile

fi

	echo -e "Adapter sequences for forward and reverse reads (tag: "${tagname}") have been created."

