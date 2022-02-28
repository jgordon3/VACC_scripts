#!/bin/sh

#  Refseq_intersects.sh
#  
#
#  Created by Jonathan Gordon on 4/21/14.
#
# REQUIRES: awk script ACC_to_GNAME.awk
# REQUIRES: kgXref table

FOLDER=/Users/jgordon3/Desktop/intersects
FILES="$FOLDER/mm10*"
TEST=NO_DAY00_TRUNC.bed
TESTNAME=$(echo $TEST | awk -F "." '{print "COUNTS_" $1 ".txt"}')
if [ -e $TESTNAME ]; then rm $TESTNAME; fi
printf "%b\t" "Genomic_region" "Intersects" "Unique_accessions" "Unique_genes" "\n" >>$TESTNAME
sort -k6n kgXref_04_22_2014 > sorted_kgXref

for f in $FILES; do
    PREFIX=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $3}')
    OUTPUT="INTERSECT_""$PREFIX""_""$TEST"
    if [[ "$PREFIX" = "introns" ]] || [[ "$PREFIX" = "codingexons" ]]; then
        intersectBed -wb -a $TEST -b $f | awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$7,$8}' | sort -k4n > $OUTPUT
    else
        intersectBed -wb -a $TEST -b $f | awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$7}' | sort -k4n  > $OUTPUT
    fi
    COUNT=$(cat $OUTPUT | wc -l);
    UNIQUEACC=$(cat $OUTPUT | awk '{ if (a[$4]++ == 0) print $0; }' "$@" | wc -l)
    awk -f ACC_to_GNAME.awk $OUTPUT > tempx.txt
if [[ "$PREFIX" = "introns" ]]; then
    UNIQUENAME=$(cat tempx.txt| awk '{ if (a[$6]++ == 0) print $0; }' "$@" | wc -l); rm tempx.txt
    printf "%b\t" "$PREFIX" "$COUNT" "$UNIQUEACC" "$UNIQUENAME" "\n" >> $TESTNAME
    awk '{ if ( $5 == "intron_0") print $0; }' $OUTPUT > temp1.txt
    FIRSTINTRONS=$(cat temp1.txt | wc -l)
    UNIQUEINTRONS=$(cat temp1.txt | awk '{ if (a[$4]++ == 0) print $0; }' "$@" | wc -l)
    printf "%b\t" "First_intron" "$FIRSTINTRONS" "$UNIQUEINTRONS" "\n" >> $TESTNAME
    rm temp1.txt
elif [[ "$PREFIX" = "codingexons" ]]; then
    UNIQUENAME=$(cat tempx.txt| awk '{ if (a[$6]++ == 0) print $0; }' "$@" | wc -l); rm tempx.txt
    printf "%b\t" "$PREFIX" "$COUNT" "$UNIQUEACC" "$UNIQUENAME" "\n" >> $TESTNAME
    awk '{ if ( $5 == "cds_0") print $0; }' $OUTPUT > temp2.txt
    FIRSTEXONS=$(cat temp2.txt | wc -l)
    UNIQUEEXONS=$(cat temp2.txt | awk '{ if (a[$4]++ == 0) print $0; }' "$@" | wc -l)
    printf "%b\t" "First_exon" "$FIRSTEXONS" "$UNIQUEEXONS"  "\n" >> $TESTNAME
    rm temp2.txt
else
    UNIQUENAME=$(cat tempx.txt| awk '{ if (a[$5]++ == 0) print $0; }' "$@" | wc -l); rm tempx.txt
    printf "%b\t" "$PREFIX" "$COUNT" "$UNIQUEACC" "$UNIQUENAME" "\n" >> $TESTNAME
    fi
done

