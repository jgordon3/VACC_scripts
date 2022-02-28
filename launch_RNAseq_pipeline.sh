#!/bin/bash

#  launch_RNAseq_pipeline.sh
#  
#
#  Created by Jonathan Gordon on 12/12/15.
#
DIVIDER=`eval printf '=%.0s' {1..100}`

usage ()
{
printf "%s\n" "This script launches entire pipeline for a whole RNAseq experiment"
printf "%s\n" "Minimum requirements are a tab-deliniated experiment file with all fastqs"
printf "%s\n" "and short names seperated by a tab. Example:"
printf "%s\n" "HWR1D14_ATGTCA_L003_R1_000.fastq.gz  MSC_D14_HWR1"
printf "%s\n" "HWR2D14_ATGTCA_L003_R1_000.fastq.gz  MSC_D14_HWR2"
printf "%s\n" "HWR3D14_ATGTCA_L003_R1_000.fastq.gz  MSC_D14_HWR3" ""
printf "%s\n" "DO NOT INCLUDE THE MATE PAIR FOR PAIRED_END IN THIS LIST" ""
printf "%s\n" "You should create an experiment folder containing all fastqs (and experiment_file.txt)"
printf "%s\n" "and submit this script from that folder, otherwise full paths are needed in"
printf "%s\n" "experiment_file.txt." ""
printf "%s\n" "Program specific options are indicated by flags"
printf "%s\n" "$DIVIDER" ""
printf "%s\n" "CUFFLINKS: run_tophat2_hisat2_on_SGE.sh" ""
printf "%s\n" "-G Genome for alignment. Default is hg38. Other options are hg19, mm9, mm10"
printf "%s\n" "-p: run is paired-end data and will search for matched read named _R2_"
printf "%s\n" "-t: trim fastq with trimmomatic prior to alignment"
printf "%s\n" "-k: Kind mode. Will limit the job submissions to only allow 3 TOPHATs at a time"
printf "%s\n" "allowing for 4 open processors"
printf "%s\n" "" "$DIVIDER" ""
exit 1;
}

while getopts "E:ptG:h" opt; do
    case $opt in
        E) EXP_FILE=$OPTARG;;          #fastq is required
        p) PAIRED="-p";;
        G) GENOME=$OPTARG;;
        t) TRIM="-t";;
        k) KIND=1;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))
if [[ $HELP = 1 ]]; then usage; fi

#EXPERIMENT FILE
if [ ! -e $EXP_FILE ]; then usage; fi
FILE_LIST=$(awk 'BEGIN {OFS=" "; FS="\t"}{print $1}' $EXP_FILE)

#Check if conditions are unique


# check KIND
# check other options
if [ ! -z $PAIRED ]; then TOP_OPTIONS="$PAIRED"; fi
if [ ! -z $TRIM ] && [ ! -z $TOP_OPTIONS ]; then TOP_OPTIONS="$TOP_OPTIONS $TRIM"; else TOP_OPTIONS="$TRIM"; fi

JOBIDS=""
for FILE in $FILE_LIST; do
    JOB=$(qsub ~/scripts/run_tophat2_hisat2_on_SGE.sh -F $FILE -G $GENOME $TOP_OPTIONS);
    JOBID=$(awk -v RS=[0-9]+ '{print RT+0;exit}' <<< "$JOB");
    JOBIDS="$JOBIDS $JOBID"
done





