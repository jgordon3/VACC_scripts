#!/bin/bash
#
# Exports UCSC files
#


#bash /slipstream/galaxy/production/galaxy-dist/tools/custom_tools/macs2_bdgcmp.sh /slipstream/galaxy/production/data/files/042/dataset_42859_files logFE /slipstream/galaxy/production/data/files/042/dataset_42868.dat

UCSC_dir=/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/

dir=$1
file_type=$2
OUT=$3

t=$(echo $dir/*treat_pileup.bdg)
c=$(echo $dir/*control_lambda.bdg)

if [ $file_type == "logFE" ] | [ $file_type == "FE" ]; then
    p=1
    else
    p=0
fi

if [ ! -f $t ] || [ ! -f $c ]; then
    echo "waiting for files to move..."
    sleep 30
fi

