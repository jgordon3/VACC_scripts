#!/bin/sh

#  Rename.sh
#  
#
#  Created by Jonathan Gordon on 9/19/13.
#
FILEPATH=/Users/jgordon3/Desktop/Corr/*.bb

for f in *; do
# let me know it is working...
echo "Processing $f file..."
FILE=$(echo "$f" | awk -F '/' '{print $NR}')
FILE=$(echo "$FILE" | awk -F '.' '{print $1}')
EXT=$(echo "$f" | awk -F '.' '{print $2}')
if [ EXT = "bowtie" ]; then echo ""; else EXT=$(echo "$f" | awk -F '.' '{print $2"."$3}'); fi
LONGDESC=$(echo "$FILE" | awk -F '_' '{print toupper ($2 "_" $1 $3 "_"$4 "_" $6)}')
ALL="$LONGDESC.$EXT"
mv $f $ALL
done
