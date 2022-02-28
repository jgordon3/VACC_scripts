#!/bin/bash
#$ -N SAMTOOLS
#$ -cwd
#$ -j y
# Request [num] amount of [type] nodes
#$ -pe threads 1

######################## DEFAULTS #############################################

DATE=`date`; NAME="J. Gordon"; MODE=1; FILE="BAM";

######################## USAGE ###############################################

usage ()
{
    printf "%s\n" "######################## USAGE ###############################################"
    printf "%s\t" "SAMTOOLS"
    printf "\n"
    printf "%s\n" "REQUIRES (NON OPTIONAL): A BAM or SAM file"
    printf "%s\n" "-F flag switches from SAM to BAM (default)"
    printf "%s\n" "OPTIONS: -M: [MODE]: 1 or V: (default) SAMTOOLS view: SAM<->BAM conversion "
    printf "%s\n" "OPTIONS: -M: [MODE]: 2 or S: SAMTOOLS sort: sort alignment file"
    printf "%s\n" "OPTIONS: -M: [MODE]: 3 or I: SAMTOOLS index: index alignment"
    printf "%s\n" "OPTIONS: -M: [MODE]: 4 or SI: SAMTOOLS sort and index"
    printf "%s\n" "example: samtools.sh -F SAM -M 3 bamfile"
exit 1;
}

######################## FLAGS ##################################################

while getopts "F:M:h:" opt; do
    case $opt in
        F) FILE=$OPTARG;;
        M) MODE=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

PATH_FILE=$(readlink -f $1)
if [ ! -f "$PATH_FILE" ]; then echo "Need valid path to bam of sam file"; fi

# write sam checks here

########################## FILE NAMES ########################################################
BAM_FILE=`basename $PATH_FILE`
ABRV_NAME=${BAM_FILE%.bam}
# Check to see if file came through STAR pipeline to shorten name
CHECK_STAR=${ABRV_NAME%_STAR_ALIGNMENT}
if [[ $ABRV_NAME == $CHECK_STAR ]]; then ABRV_NAME=$ABRV_NAME; else ABRV_NAME=$CHECK_STAR; fi
SORT_OUT=$ABRV_NAME"_SORTED_ALIGNMENT"
INDEX_OUT=$ABRV_NAME"_INDEX.bai"

if [ "$MODE" = "1" ] | [ "$MODE" = "V" ]; then
    samtools view
elif [ "$MODE" = "2" ] | [ "$MODE" = "S" ]; then
    samtools sort $BAM_FILE $SORT_OUT
    rm $BAM_FILE
elif [ "$MODE" = "3" ] | [ "$MODE" = "I" ]; then
    samtools index $BAM_FILE $INDEX_OUT
elif [ "$MODE" = "4" ] | [ "$MODE" = "SI" ]; then
    samtools sort $BAM_FILE $SORT_OUT
    rm $BAM_FILE
    INDEX_OUT=$SORT_OUT"_INDEX.bai"
    SORT_OUT=$SORT_OUT".bam"
    samtools index $SORT_OUT $INDEX_OUT
fi


###### TROUBLESHOOT
#echo "MODE: $MODE"; echo "PATH_FILE; $PATH_FILE"; echo "ABRV_NAME; $ABRV_NAME"; echo "SORT_OUT: $SORT_OUT"


