#!/bin/sh

#  submit_RNAseq_pipeline.sh
# SUBMITS TOPHAT-CUFFLINKS-CUFFDIFF and others

while getopts "w:f:hs::r::l:c" opt; do
    case $opt in
        w) WD=$OPTARG;;
        s) SPECIES=$OPTARG;;
        r) REFERENCE=$OPTARG;;
        f) FASTQ=$OPTARG;;
        l) LIBTYPE=$OPTARG;;
        c) FNC="Y";;
        h) HELP=1;;
        \?) HELP=1;;
        *) HELP=1;;
    esac
done

shift $(($OPTIND -1))

if [ -z $WD ]; then printf "\nNeed a valid working directory (-w)"; HELP=1; else if [ -d "$WD" ]; then cd $WD; else mkdir $WD; cd $WD; fi; fi;
if [ -z $FASTQ ]; then printf "\nNeed a (full) path to fastq files (-f)"; HELP=1; else FASTQS=${FASTQ%/}; FASTQS="$FASTQS/*fastq"; fi
if ls $FASTQS &> /dev/null; then echo "fastq files are there"; else echo "Folder $FASTQ does not appear to have fastq files"; HELP=1; fi

if [ -z $HELP ]; then echo ""; else
printf "\n\nHELP FILE\n"
printf "You are probably reading this because you did something wrong."
printf "\n\n\n"
printf "This script submits several scripts to run the entire RNAseq pipeline.\n" "REQUIRED: (-w) A working directory. If the working directory does not exist it will be created\n"
printf "REQUIRED (-f) Folder containing fastq files. These can be gzipped\n"
printf 'REQUIRED: (-s) Species to align. Options are "mm10" or "hg19" The default is mm10 if unflagged'"\n"
printf 'REQUIRED: (-r) Reference to align. Options are "UCSC" or "NCBI" The default is UCSC if unflagged'
printf "\n"
printf 'OPTIONS: (-l) Library type. Options are: "fr-unstranded" (default if unflagged) "fr-firststrand" or "fr-secondstrand"'
printf "\n"
printf '         See: Trapnell, C et al (2013). Differential analysis of gene regulation at transcript resolution with RNA-seq.'
printf "\n"
printf '          Nature Biotechnology, 31(1), 46–53 for more detail on these data types'
printf "\n"
printf 'OPTIONS: (-c) File naming convention: Fastq files are in the format: "CELL_CONDITION_REPLICATE"'
printf "\n"
printf '          where at least one of these is unique. An example would be "BMSC_DAY00_R1"'
printf "\n"
printf "OPTIONS; (-h) This information"
printf "\n";
printf "This script sleeps for a while so you may want to use nohup"
exit 1; fi

SCRIPT=
# use a tempfile to count how many jobs are submitted
TEMPFILE=~/.tmp; echo 0 > $TEMPFILE;

for f in $FASTQS; do

#LOOP COUNTER
COUNT=$[$(cat $TEMPFILE) + 1]; echo $COUNT > $TEMPFILE
NAME="TOPHAT_PIPELINE""$COUNT"
echo "qsub $NAME"
#qsub -v PATHTOFASTQ=$f,SPECIES=$SPECIES,REFERENCE=$REFERENCE,FNC=$FNC,WD=$WD,LIBTYPE=$LIBTYPE -N $NAME $SCRIPT
done
COUNT=$(cat $TEMPFILE)
echo "waiting for $COUNT files"


#if [[ -f .ssh/id_rsa && -f .ssh/id_rsa.pub ]]; then check to see if both files exist


#sleep 60
#check that all the jobs completed

#LOOPCOUNT=1
#while [ $LOOPCOUNT -le $COUNT ]
#do
#printf "loop iteration: $LOOPCOUNT\n"

#LOOPCOUNT=$[$LOOPCOUNT+1]
#done
#echo "EOF"

unlink $TEMPFILE

