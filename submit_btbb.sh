#!/bin/sh

#general submit script
#OPTIONS

SPECIES="hg19"
SCRIPT="/users/j/a/jargordo/bin/scripts/run_bedToBigbed_on_VACC.sh"
FOLDER="/users/j/a/jargordo/scratch/TM_R2_PEAK_CALL_MACS2_03_19_2014/*.bed"
WD="/users/j/a/jargordo/scratch/TM_R2_PEAK_CALL_MACS2_03_19_2014"
#input files to submit

#whole folder

#for f in $FOLDER; do
for f in $FOLDER; do
qsub -v INPUT=$f,SPECIES=$SPECIES,WD=$WD $SCRIPT
done
