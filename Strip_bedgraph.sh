#!/bin/sh

#  Strip_bedgraph.sh
#  
#
#  Created by Jonathan Gordon on 4/23/14.
#
PEAKFILE=
BEDGRAPH=
file="DAY00_v_DAY07_v_DAY14_v_DAY21_INTERSECT_OL50.bed"
while read line
do
chr=$(echo $line | awk -F " " '{print $1}')
start=$(echo $line | awk -F " " '{print $2}')
stop=$(echo $line | awk -F " " '{print $3}')
peak=$(echo $line | awk -F " " '{print $4}')
PREVCHR="chr1"
if [[ "$chr" = "$PREVCHR" ]]; then
    printf "."
else
    printf "$chr"
fi
PREVCHR=$chr
awk -F "\t" -v chr="$chr" -v start="$start" -v stop="$stop" -v peak="$peak" 'BEGIN{OFS="\t";} { if ($1 == chr && $2 >= start  && $3 <= stop) print peak,$1,$2,$3,$4}' BMSC_RUNX2D00_COMBREPS_42235272_treat_pileup.bdg >> newtest.out
done <"$file"