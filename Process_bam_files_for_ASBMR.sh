#!/bin/bash

#  Process_bam_files_for_ASBMR.sh
#  
#
#  Created by Jonathan Gordon on 8/22/16.
#
EXPNAME=$1
BAM=$EXPNAME"_Aligned.out.bam"


outbam=$name"_scaled.bam"
sortedbam=$name"_scaled_sorted.bam"
fastq1=$name"_scaled_P1.fastq"
fastq2=$name"_scaled_P2.fastq"

samtools index $BAM

for f in {1..5}; do
    samtools view -b $BAM chr$f > in_chr$f.bam
done

samtools merge -f $outbam in_chr*.bam
rm in_chr*.bam
#samtools sort -n $outbam $sortedbam
bedtools bamtofastq -i $sortedbam -fq1 $fastq1 -fq2 $fastq2
cp $fastq1 $fastq2 $parent
cd $parent
done