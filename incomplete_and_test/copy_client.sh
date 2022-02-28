#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=1gb,pvmem=1gb
#PBS -l walltime=30:00:00
#PBS -N CopyClient
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

#USER INPUT
#change this to the full path to the Lane designations:
PATHTOFILES="/users/j/a/jargordo/scratch/users/r/d/rdevins/scratch/MPS/130613_SNL128_0059_AC2668ACXX/Unaligned_Lane2"
COPYPATH=

#
cd $PATHTOFILES
PROJECTFOLDER=$(ls $WD | grep Project)
cd $PROJECTFOLDER
SAMPLES="*"
echo $SAMPLES

for s in $SAMPLES; do
cd "$s"
FILES="*"
echo "$FILES"

for f in $FILES; do
FILETYPE=$(echo "$f"| awk -F '.' '{print $2 "." $3}')
if [ $FILETYPE == 'fastq.gz' ]
then
FILENAME=$(echo "$f"| awk -F '.' '{print $1}')
echo "$FILENAME is a gzip"
READ=$(echo "$FILENAME"| awk -F '_' '{print $(NF-1)}')
if [ $READ == 'R1' ] #check read designation (R1 is forward read, R2 is paired end)
then
echo "$f is replicate one"
        #cat $f >> "$COPYPATH$FILENAME.fastq.gz"
else
echo "should be R2 here: $f"
        #cat $f >> "$COPYPATH$FILENAME.fastq.gz"
fi
else
#other file
echo "$f is not gzip"
fi
done

cd "$PATHTOFILES/$PROJECTFOLDER"
done
