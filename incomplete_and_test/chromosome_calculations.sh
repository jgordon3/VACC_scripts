#!/bin/sh

#  chromosome_calculations.sh
#  
#
#  Created by Jonathan Gordon on 8/27/13.
#

for sample in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY; do
awk -v r=$sample '{if ($1==r) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7;}' gene_lengths_by_chromosome_fixed.bed > temp.1

echo $sample >> chrome_output
awk '{ sum += $7 } END { print sum }' temp.1 >> chrome_output
done