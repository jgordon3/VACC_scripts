#!/bin/sh

#  pair.sh
#  
#
#  Created by Jonathan Gordon on 6/10/16.
#
for i1 in *_R1_000.fastq.gz; do
i2=${i1/R1_000/R2_000}
paired=${i1/R1_000/paired}
echo "$i1 $i2"
done

for i1 in *_R1_000.fastq.gz; do
i2=${i1/R1_000/R2_000}
echo "$i1 $i2"
done

for i1 in *_R1_000.fastq.gz; do
i2=${i1/R1_000/R2_000}
qsub run_kallisto_test.sh $i1 $i2
done