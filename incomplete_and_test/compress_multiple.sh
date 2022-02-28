#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash


#  compress_multiple.sh
#  
#
#  Created by Jonathan Gordon on 11/20/15.
#

file=$1


for f in $file*; do
    x=${f%_sequence.txt};
    y=$x"_GAII_sequence.fastq";
    #mv $f $y
    #gzip $y
    z=$y".gz"

done


