#!/bin/sh


PARENTDIR=/users/j/a/jargordo/scratch/TOPHAT_03_19_2014
FILES=$PARENTDIR/*
CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/mm10.chrom.sizes"

for f in $FILES; do

qsub -v PARENTDIR=$PARENTDIR,FOLDER=$f,CHROMSIZE=$CHROMSIZE ~/bin/scripts/run_bam_to_bigwig.sh
done