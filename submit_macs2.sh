#!/bin/sh

#submit script for run_MACS2_on_VACC.sh
#OPTIONS
GENOME="hg19"
#GENOME="mm10"

# SET IDR PARAMETER
# IDR=0 means MACS will call just default PVALUE; IDR=1 will call just the IDR peaks (1e-03); IDR=2 will call both sets of peaks (set Pvalue, 1e-03)
IDR=2
# SET PVALUE AND MFOLD
PVALUE=1e-10
MFOLDL="5"
MFOLDU="50"
#SET TO 0 BEGRAPH WILL BE MADE WITH --SPMR, IF SET TO 1 BEDGRAPH WILL NOT BE MADE (FOR IDR ANALYSIS)
NOSPMR=0
WD=/users/j/a/jargordo/scratch/TM_MCF_MDA_HISTONE_TREATMENTS_MACS_06_06_2014

# SET BOWTIE FOLDERS
BOWTIEFOLDER=/users/j/a/jargordo/scratch/TM_MCF_MDA_HISTONE_TREATMENTS_BOWTIE_05_09_2014
# IF YOU WANT TO DO REPLICATE ANALYSIS
BOWTIE2FOLDER=



BOWTIEFILES="$BOWTIEFOLDER/*bowtie"
for f in $BOWTIEFILES; do

FILE=$(echo "$f" | awk -F '/' '{print $NF}')
NAME=$(echo "$FILE" | awk -F '_' '{print $1 "_" $2 "_R2"}')
BOWTIE2=$BOWTIE2FOLDER/"$NAME""*.bowtie"

echo $f

qsub -v BOWTIE=$f,GENOME=$GENOME,WD=$WD,IDR=$IDR,PVALUE=$PVALUE,MFOLDU=$MFOLDU,MFOLDL=$MFOLDL,NOSPMR=$NOSPMR /users/j/a/jargordo/bin/scripts/run_MACS2_on_VACC.sh
done
