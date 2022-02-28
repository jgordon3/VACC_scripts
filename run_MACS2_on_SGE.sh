#!/bin/bash
#$ -N MACS2
#$ -cwd
#$ -j y
#######################################################
####################### DEFAULTS GO HERE ##############
#######################################################

BOWTIE_VERSION=`macs2 --version`; COMPRESS="y"; DATE=`date`; FILTER="y"; GENOME="hg38";
NAME="J. Gordon"; RECIPE="SE"; PVALUE="1e-05"; SAVEBED="y"; ADD_PEAKCALL="n"; MFOLD="5 50"; MACSFOLDER=`pwd`;

####### HELP #########

usage ()
{
printf "%s\n" "run_MACS2_on_SGE.sh" ""
printf "%s\n" "REQUIRES (NON OPTIONAL): A treatment and input BAM file."
printf "%s\n" "This can be providied by the -T flag with a path/to/a/treatment.bam and -I flag with a path/to/a/control.bam "
printf "%s\n" "REQUIRES: A reference genome"
printf "%s\n" "The script will default to hg38 or the -g flag can be used. Options are hg38, hg19, mm9, mm10" ""
printf "%s\n" "OPTIONS:" "-p: MACS p value"
printf "%s\n" "OPTIONS:" "-m: M Fold for MACS2"
exit 1;
}


# OPTIONS
while getopts "T:I:g:p:m:u:s:a:h" opt; do
    case $opt in
        T) TREAT=$OPTARG;;
        I) INPUT=$OPTARG;;
        g) GENOME=$OPTARG;;
        p) PVALUE=$OPTARG;;
        m) MFOLD=$OPTARG;;
        u) NAME=$OPTARG;;
        s) SAVEBED=$OPTARG;;
        a) ADD_PEAKCALL=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

#CHECK PATHS
TREAT_PATH=$(readlink -f $TREAT)
INDEX_PATH=$(readlink -f $INDEX)
if [[ ! -e $TREAT_PATH ]]; then echo "Need valid path to BOWTIE BAM/SAM file"; usage; fi
if [[ ! -e $INPUT_PATH ]]; then echo "No Input set"; fi

# CHECK FOR GZIP
EXT=${TREAT_PATH##*.}
if [ $EXT = "gz" ]; then gunzip $TREAT_PATH; TREAT_PATH=${TREAT_PATH%.gz}; fi
EXT=${INPUT_PATH##*.}
if [ $EXT = "gz" ]; then gunzip $INDEX_PATH; INDEX_PATH=${INDEX_PATH%.gz}; fi

# CHECK AND MAKE FOLDERS
TREATABRVNAME=${TREAT%.bam}
INPUTABRVNAME=${INPUT%.bam}
OUTPUT_DIR=./MACS/$TREATABRVNAME
if [ ! -d "./MACS" ]; then mkdir MACS; fi
if [ ! -d "$OUTPUT_DIR" ]; then mkdir $OUTPUT_DIR; fi
cd $OUTPUT_DIR

# SET GENOME FILES
if [ "$GENOME" = "hg38" ] || [ "$GENOME" = "hg19" ]; then GEN="hs";
    elif [ "$GENOME" = "mm9" ] || [ "$GENOME" = "mm10" ]; then GEN="mm";
fi  #hs 2.7e9 mm 1.87e9

# get chromosome sizes
CHROMSIZE=$GENOME.chrom.sizes
fetchChromSizes $GENOME > $CHROMSIZE

###################   RUN MACS2 HERE ##############################

macs2 callpeak -B -t $TREAT_PATH -c $INPUT_PATH -m $MFOLD -f BAM -n $TREATABRVNAME -g $GEN -p $PVALUE

# run MACS2 for fold enrichment
macs2 bdgcmp -t "${TREATABRVNAME}_treat_pileup.bdg" -c "${TREATABRVNAME}_control_lambda.bdg" -o "${TREATABRVNAME}_FE.bdg" -m FE

# run MACS2 for logLR
macs2 bdgcmp -t "${TREATABRVNAME}_treat_pileup.bdg" -c "${TREATABRVNAME}_control_lambda.bdg" -o "${TREATABRVNAME}_logLR.bdg" -m logLR -p 0.00001

# Convert bedgraphs to bigwig
wigToBigWig -clip "${TREATABRVNAME}_treat_pileup.bdg" $CHROMSIZE "${TREATABRVNAME}_treat_pileup.bw"

#run others

#Covert beds to bigbed

# Calculate coverage
MAPPED_READS=`samtools view -c -F 260 $TREAT_PATH`
READS_IN_PEAKS=`bedtools coverage -counts -abam $TREAT_PATH -b "${TREATABRVNAME}_peaks.narrowPeak" | awk -F"\t" '{sum += $11} END {print sum}'`
FRAC_IN_PEAKS=`echo "scale=4 ; $READS_IN_PEAKS / $MAPPED_READS" | bc`
NUM_OF_PEAKS=`wc -l ${TREATABRVNAME}_peaks.narrowPeak`

if [ ! -f "../MACS_summary.txt" ]; then 

