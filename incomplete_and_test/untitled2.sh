#!/bin/sh

# We find all files in path, feed them into ls -l with xargs,
# and sort them on the size column.
# We can't depend on ls' own sort when using xargs since enough
# files will end up splitting between several ls calls.
# Then we read the lines in order, and check for duplicate sizes.

FILE=*
echo $FILE | xargs ls -l | sort -k 5,6 |
while read PERMS LINKS USER GROUP SIZE M D Y FILE
do
    if [ "$SIZE" -eq "$LASTSIZE" ]; then
        echo "$LASTFILE $FILE"; echo "$LASTFILE $FILE" >> pairs.txt
    else
        LASTSIZE="$SIZE"; LASTFILE="$FILE"
    fi
done
