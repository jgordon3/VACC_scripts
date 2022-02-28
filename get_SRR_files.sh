#!/bin/bash

#  get_SRR_files.sh
#  
#
#  Created by Jonathan Gordon on 12/15/15.
#



CONFIG=$1
perl -pi -e 's/\r/\n/g' $CONFIG #fix MAC tsv

SRR=($(cut -f 1 $CONFIG))
GSM=($(cut -f 2 $CONFIG))
LONG_NAME=($(cut -f 3 $CONFIG))
i=0
len=${#SRR[@]}

while [ $i -lt $len ]
    do #echo $i
        srr=${SRR[$i]}
        name=${LONG_NAME[$i]}
        name=${name/_/}; name=${name/_/}; 
        name=$(echo $name | awk '{print toupper($0)}')"_"
        gsm=${GSM[$i]}
        echo $srr $name $gsm
        filename="$name"_"$gsm".fastq
        echo $filename
        i=$(($i + 1))

done
