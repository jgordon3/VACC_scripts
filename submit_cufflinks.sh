#!/bin/sh

#general submit script
#OPTIONS

SPECIES="mm10"
SCRIPT="/users/j/a/jargordo/bin/scripts/run_CUFFLINKS_on_VACC.sh"
FOLDER="/users/j/a/jargordo/scratch/CT_TOPHAT_04_15_2014"
FILES="$FOLDER/*"
WD="/users/j/a/jargordo/scratch/CT_TOPHAT_04_15_2014"
REFERENCE="UCSC"


for f in $FILES; do
#SHOULD BE A DIRECTORY FROM TOPHAT
if [ -d $f ]; then
    if [[ $f = "$FOLDER/BIGWIGS" ]]; then echo "";
    else qsub -v INPUT=$f,SPECIES=$SPECIES,WD=$WD,REFERENCE=$REFERENCE $SCRIPT;
    fi
fi
done
