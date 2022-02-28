#!/bin/bash

# breaking down the awk scripts from the Nature methods paper
#awk '(NR-2)%4==0' s_4_sequence.txt > out.txt #returns every 4th line starting at line 1

#reads starting at line one increments total and "read"
cat s_4_sequence.txt | awk '((NR-2)%4==0){read=$1;total++;count[read]++}END{for(read in count){if(count[read]==1){unique++}} print "total reads: "total, "unique reads: "unique, "percent unique: "unique*100/total}'

#clasic min max awk
#awk ‘{if(min==""){min=max=$1}; if($1>max) {max=$1}; if($1< min) {min=$1}; total+=$1; count+=1} END {print total/count, min, max}’ data.txt