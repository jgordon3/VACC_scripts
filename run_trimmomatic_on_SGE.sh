#!/bin/sh
#  run_trimmomatic_on_SGE.sh
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 4

############################# DEFAULTS ######################
TRIMMOMATIC_DIR=/slipstream/home/jonathan/bin/Trimmomatic-0.35
TRIMMOMATIC_ADAPTERS_DIR=/slipstream/home/jonathan/bin/Trimmomatic-0.35/adapters/TruSeq3-SE.fa
RECIPE="SE"; #-R|--recipe
ILLUMINACLIP_FILE=/slipstream/home/jonathan/bin/Trimmomatic-0.35/adapters/TruSeq3-SE.fa: #-i|--illumina_clip_file
CLIP_LEADING=3 #-l|--clip_leading
CLIP_TRAILING=3 #-t|--clip_trailing
WINDOW="4:20" #-w|--window
MINLEN=36 #-m | --min_length

############################# FLAGS #########################

while getopts "F:G:t:u:i:l:w:R:h:" opt; do
    case $opt in
        F) FASTQ=$OPTARG;;
        i) ILLUMINACLIP_FILE=$OPTARG;;
        l) CLIP_LEADING=$OPTARG;;
        t) CLIP_TRAILING=$OPTARG;;
        w) WINDOW=$OPTARG;;
        m) MINLEN=$OPTARG;;
        R) RECIPE=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

ABRV_NAME=${FASTQ%.fastq}
FILT_FASTQ=$ABRV_NAME"_trim_filt.fastq"

java -Djava.io.tmpdir=/slipstream/home/jonathan/bin/tmp -jar $TRIMMOMATIC_DIR/trimmomatic-0.35.jar SE -threads 4 -phred33 $FASTQ $FILT_FASTQ ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS_DIR:2:30:10 LEADING:$CLIP_LEADING TRAILING:$CLIP_TRAILING SLIDINGWINDOW:$WINDOW MINLEN:$MINLEN;