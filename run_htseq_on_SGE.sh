#!/bin/bash
#$ -N HTSEQ
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 1

# don't change this unless absolutely needed
DATE=`date`

while getopts "D:BG:h" opt; do
    case $opt in
        D) TOPDIR=$OPTARG;; #fastq is required
        B) BAM_FILE=$OPTARG;;
        G) GTF_FILE=$OPTARG;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))


cd $TOPDIR
if [ ! -z $BAM_FILE ]; then printf "%s\n" "$BAM_FILE"; else BAM_FILE=$(find . -name *accepted* -print); BAM_FILE=$(basename $BAM_FILE); echo $BAM_FILE; fi;


ABRV_NAME=${TOPDIR%_TOPHAT_OUT}; echo $ABRV_NAME
OUTPUT_NAME=$ABRV_NAME"_htseq_counts.txt";

htseq-count -f bam $BAM_FILE $GTF_FILE > $OUTPUT_NAME