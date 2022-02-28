#!/bin/bash
#$ -N WAIT
#$ -cwd
#$ -j y

STARTTIME=$(date +%s)

for i in $(seq 1 10); do
    sleep 120
    CURRENTTIME=$(date +%s)
    ELAPSED=$(($CURRENTTIME - $STARTTIME))
    echo "Waiting for $FILE for $ELAPSED seconds"
done
