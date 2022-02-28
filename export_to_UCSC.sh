#!/bin/sh
#export all data from MACS2 html output

dir=$1 #MACS2 tag directory

#4 data files

a=$(echo $dir/*control_lambda.bdg)
b=$(echo $dir/*peaks.narrowPeak)
c=$(echo $dir/*summits.bed)
d=$(echo $dir/*treat_pileup.bdg)

counter=0
while [$counter -lt 2 ]; do
    if [ ! -f $a ] || [ ! -f $b ] || [ ! -f $c ] || [ ! -f $d ]; then
        echo "waiting for files to move..."
        sleep 20
    fi
counter=$counter+1
done


cmd="macs2 bdgcmp -t $t -c $c -m $met -o $OUT -p $p"
echo $cmd
$cmd 2>&1
