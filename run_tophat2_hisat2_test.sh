#!/bin/bash
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 4

#GENERAL STUFF
DIVIDER=`eval printf '=%.0s' {1..100}`
TOPHAT_VERSION=`tophat2 --version`
DATE=`date`

usage ()
{
case $ERROR in
    1) MESSAGE="You did not supply a fastq with the -F flag or it does not exist";;
    2) MESSAGE="Paired end flagged (-p) but no mate file";;
    3) MESSAGE="Genome was not set ... which is weird because it is supposed to default to hg38. Please set it with -G";;
    *) MESSAGE="Help file activated or Random Error: resubmit";;
esac
printf "%s\n" "$DIVIDER" "" "$MESSAGE" "" "This is script is now submitted via qsub with the command" "" "qsub run_tophat2_hisat2_on_SGE.sh -F /path/to/fastq" "" "****  optional flags: ******** " "" "-p: Paired end data and will search for matched read (R2) unless a path to the mate is supplied. Default assumes no pair" "-g: Genome for alignment. Default is hg38. Some other options are hg19, mm9, mm10 ... etc" "-t: [y/n]: filter fastq with Trimmomatic for Illumina Adapters, quality, and length. Default is y" "-u: User info for README file" "-q: [y/n] Submit fastqc script for raw and/or trimmed fastq files. Default is n." "$DIVIDER" ""
exit 1;
}

#set flags here

while getopts "F:pG:h" opt; do
    case $opt in
        F) FASTQ=$OPTARG;;          #fastq is required
        p) PAIREDEND=1;;
        G) GENOME=$OPTARG;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))

printf "%s\n" "$DIVIDER" ""

if [ -e $FASTQ ]; then printf "%s\n" "Fastq: $FASTQ"; else ERROR=1; usage; fi
if [[ $HELP = 1 ]]; then usage; fi
if [ -z $GENOME ]; then GENOME="mm10"; printf "%s\n" "Genome set by default to: $GENOME"; else GENOME=$GENOME; printf "%s\n" "Genome set by flag (-G) to: $GENOME"; fi

#CHECK FOR PAIRED FILE
if [[ $PAIREDEND = 1 ]]; then FASTQ2=${FASTQ//_R1_/_R2_}; TYPE="PE100"; else TYPE="SE100"; fi  #fix this to be more general # check to see if there is more than one R1 in file name
if [[ $PAIREDEND = 1 ]] && [ ! -e $FASTQ2 ]; then ERROR=2; usage; else printf "%s\n" "Paired fastq was found: $FASTQ2"; fi

#CHECK FOR GZIP
EXT=${FASTQ##*.}
if [ "$EXT" = "gz" ]; then gunzip $FASTQ; gunzip $FASTQ2; FASTQ=${FASTQ%.gz}; FASTQ2=${FASTQ2%.gz}; printf "%s\n" "Fastq is gzipped ..... unzipping"; else printf "%s\n" "Not gzipped.....doing nothing"; fi

#TRIM

#OUPUT
ABRV_NAME=${FASTQ%.fastq}
OUTPUT_NAME=$ABRV_NAME"_TOPHAT_OUT"; mkdir $OUTPUT_NAME;

# SET GENOME INDEX AND GTF
if [ "$GENOME" = "hg38" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/hg38/hg38canon/bowtie2_index/hg38canon;
    elif [ "$GENOME" = "hg19" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/hg19canon/bowtie2_index_canon/hg19canon;
    elif [ "$GENOME" = "mm10" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/mm10/bowtie2_index_canon/mm10canon;
    elif [ "$GENOME" = "mm9" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/mm9/bowtie2_index_canon/mm9canon;
    else ERROR=3; usage;
fi



