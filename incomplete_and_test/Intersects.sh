#!/bin/sh

#  Intersects.sh
#  
#
#  Created by Jonathan Gordon on 4/21/14.
#

FOLDER=/Users/jgordon3/Desktop/intersects
FILES="$FOLDER/*peaks_IDR_filtered.narrowPeak"
rm *_INTERSECT*
rm counts.txt

for f in $FILES; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    for s in $FILES; do
    DAY2=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    if [[ "$DAY1" = "$DAY2" ]]; then echo "same"; else
    PREFIX=$DAY1"_v_"$DAY2
    OUTPUT1="$PREFIX""_INTERSECT_UNIQUE.bed"
    intersectBed -v -a $f -b $s | awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$7}' > $OUTPUT1
    fi
    done
done

for f in $FILES; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    for s in $FILES; do
    DAY2=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    if [[ "$DAY1" = "$DAY2" ]]; then echo "same"; else
    PREFIX=$DAY1"_v_"$DAY2
    OUTPUT2="$PREFIX""_INTERSECT_OL50.bed"
    intersectBed -wb -a $f -b $s -f 0.5 |  awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$7,$20}' > $OUTPUT2
    fi
    done
done


#ROUND 2
NEWOLS="$FOLDER/*_INTERSECT_OL50.bed"
for f in $NEWOLS; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
DAY2=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $3}')
    for s in $FILES; do
    DAYX=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    if [[ "$DAYX" = "$DAY1" ]] || [[ "$DAYX" = "$DAY2" ]]; then echo ""; else
    PREFIX=$DAY1"_v_"$DAY2"_v_"$DAYX
    OUTPUT2=$PREFIX"_INTERSECT_OL50.bed"
    intersectBed -wb -a $f -b $s -f 0.5 | awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$4,$5,$12}' > $OUTPUT2
    fi
    done
done

NEWUNIQS="$FOLDER/*_INTERSECT_UNIQUE.bed"
for f in $NEWUNIQS; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
DAY2=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $3}')
    for s in $FILES; do
    DAYX=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
    if [[ "$DAYX" = "$DAY1" ]] || [[ "$DAYX" = "$DAY2" ]]; then echo ""; else
    PREFIX=$DAY1"_v_"$DAY2"_v_"$DAYX
    OUTPUT1=$PREFIX"_INTERSECT_UNIQUE.bed"
    intersectBed -v -a $f -b $s > $OUTPUT1
    fi
    done
done

#ROUND3
FINALOLS="$FOLDER/*_INTERSECT_OL50.bed"
for f in $FINALOLS; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
DAY2=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $3}')
DAY3=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $5}')
TYPE=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $7}')
echo "$DAY1 $DAY2 $DAY3 $TYPE"
if [[ "$DAY3" = DAY* ]] || [[ "$TYPE" = "OL50.bed" ]]; then
    for s in $FILES; do
        DAYX=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
        echo "DAYX= $DAYX"
        if [[ "$DAYX" = "$DAY1" ]] || [[ "$DAYX" = "$DAY2" ]] || [[ "$DAYX" = "$DAY3" ]]; then continue; else
        PREFIX1=$DAY1"_v_"$DAY2"_v_"$DAY3
        PREFIX=$PREFIX1"_v_"$DAYX
        OUTPUT1=$PREFIX"_INTERSECT_OL50.bed"
        intersectBed -wb -a $f -b $s -f 0.5 | awk 'BEGIN{OFS="\t";} {print $1,$2,$3,$4,$5,$6,$13}'  > $OUTPUT1
        fi
    done
else continue;
fi
done

FINALUNIQS="$FOLDER/*_INTERSECT_UNIQUE.bed"
for f in $FINALUNIQS; do
DAY1=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
DAY2=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $3}')
DAY3=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $5}')
TYPE=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $7}')
echo "$DAY1 $DAY2 $DAY3 $TYPE"
if [[ "$DAY3" = DAY* ]] || [[ "$TYPE" = "UNIQUE.bed" ]]; then
    for s in $FILES; do
        DAYX=$(echo $s | awk -F "/" '{print $NF}' | awk -F "_" '{print $1}')
        echo "DAYX= $DAYX"
        if [[ "$DAYX" = "$DAY1" ]] || [[ "$DAYX" = "$DAY2" ]] || [[ "$DAYX" = "$DAY3" ]]; then continue; else
        PREFIX1=$DAY1"_v_"$DAY2"_v_"$DAY3
        PREFIX=$PREFIX1"_v_"$DAYX
        OUTPUT1=$PREFIX"_INTERSECT_UNIQUE.bed"
        intersectBed -v -a $f -b $s > $OUTPUT1
        fi
    done
else continue;
fi
done


for f in *; do
COUNT=$(wc -l $f);
SPLIT=$(echo $COUNT | awk -F " " '{print $2 "\t" $1 }')
printf "$SPLIT" >>counts.txt
printf "\n" >>counts.txt
done


#FILES="mm10*"
#for f in $FILES; do
#PREFIX=$(echo $f | awk -F "/" '{print $NF}' | awk -F "_" '{print $1 "_" $2 "_" $3 "_" $4}')
#echo $PREFIX
#OUTPUT="$PREFIX""_DAY00_v_DAY07_v_DAY14_v_DAY21_INTERSECT_UNIQUE.bed"
#intersectBed -wb -a DAY00_v_DAY07_v_DAY14_v_DAY21_INTERSECT_UNIQUE.bed -b $f -f 0.5 > $OUTPUT
#done

#for f in "*INTERSECT_*"; do
#COUNT=$(wc -l $f);
#printf $COUNT >> counts.txt
#printf "\n" >> counts.txt
#done

