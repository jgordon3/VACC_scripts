#!/bin/bash
find . -type f -print0 | (
    while read -d "" FILE ; do FILES=("${FILES[@]}" "$FILE") ; done

    ls -la "${FILES[@]}" | awk '{$1=$2=$3=$4=$6=$7="";print}' > /Users/Jgordon3/Desktop/Listed_Files.txt
    ls -la "${FILES[@]}" | awk '{print $5}' | sort -k1,1nr | uniq -d > /Users/Jgordon3/Desktop/Repeated_Sizes.txt

)

awk 'BEGIN{print "Size (bytes)  Files"}FNR==NR{a[$1];next} $1 in a' Repeated_Sizes.txt Listed_Files.txt > Duplicates_Files.txt

#rm /Users/Jgordon3/Desktop/Listed_Files.txt
#rm /Users/Jgordon3/Desktop/Repeated_Sizes.txt