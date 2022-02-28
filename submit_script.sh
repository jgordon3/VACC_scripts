#!/bin/sh

#general submit script
#OPTIONS

SPECIES="mm10"
SCRIPT="/users/j/a/jargordo/bin/scripts/run_CUFFDIFF_on_VACC.sh"
WD="/users/j/a/jargordo/scratch/CT_TOPHAT_SS_03_27_2014"
REFERENCE="UCSC"
FNC="Y"

qsub -v SPECIES=$SPECIES,WD=$WD,REFERENCE=$REFERENCE,IFN=$IFN $SCRIPT

