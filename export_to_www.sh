#!/bin/bash
#Export to file to display on UCSC browser
# needs export_to_www.xml

d=$1
met=$2
OUT=$3
t=$(echo $d/*treat_pileup.bdg)
c=$(echo $d/*control_lambda.bdg)
if [ $met == "logFE" ] | [ $met == "FE" ]; then
p=1
else
p=0
fi
if [ ! -f $t ] || [ ! -f $c ]; then
echo waiting for files to move...
sleep 10
fi
if [ ! -f $t ] || [ ! -f $c ]; then
echo waiting some more...
sleep 10
fi
if [ ! -f $t ] || [ ! -f $c ]; then
echo giving up! cannot find bdg files.
exit 1
fi



cmd="macs2 bdgcmp -t $t -c $c -m $met -o $OUT -p $p"
echo $cmd
$cmd 2>&1
