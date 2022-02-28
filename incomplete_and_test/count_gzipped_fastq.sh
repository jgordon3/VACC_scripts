#!/bin/bash
# count reads in a gzipped fastq
for f in *.fastq.gz; do 
	f1=${f%.gz}; 
	gunzip $f -c >> $f1; 
	count=`wc -l $f1`; 
	fsize=`du -h $f1`; 
	lcount=${count% *}; 
	rcount=$(expr $lcount / 4); 
	printf "$f1\t$rcount\t$fsize\n" >> readcounts.txt
	rm $f1
done