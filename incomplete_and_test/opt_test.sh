#!/bin/sh

#  opt_test.sh
#  
#
#  Created by Jonathan Gordon on 3/25/14.
#A (no colon) tells that the option is not reqired.
#A : (colon character) tells that the option has a required argument.
#A :: (two consequent colon character) tells that the option has an optional argument. (if "on" requires a value)
#FLAG HANDLING
while getopts "n:vhc::f::u::" opt; do
    case $opt in
        n) NAME=$OPTARG;;
        v) VAL=1;;
        h) HELP=1;;
        c) CSEP=$OPTARG;;
        f) FILEPATH=$OPTARG;;
        u) URL=$OPTARG;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
if [ -z $NAME ]; then echo "Need a name for the control file (-n)"; exit 1; else echo "$NAME"; fi
if [[ ! -z $VAL ]]; then echo "Validate files"; fi
if [ -z $HELP ]; then echo ""; else
echo "GENERATES A UCSC GENOME BROWSER CONTROL FILE FROM EVERY FILE IN A DIRECTORY
echo 'Looks for files labeled .bb (BigBed) or .bw (Bigwig) and generates an entry';
echo 'USAGE Generate_control_file.sh -n [NAME] -f [FILEPATH] -v -c [CELL/CONDITION/REPLICATE] -u [URL]'
echo 'REQUIRES: A control file name (-n "FILENAME").';
echo 'The program will remove this file if it already exists.';
echo 'OPTION: Validate files before writing (-v).'
echo 'Uses UCSC validateFiles tool to check if files are valid format before writing a entry';
echo 'OPTION: Colour by variable (-c). If file is in FNC format ("CELL_CONDITION_REPLICATE"),';
echo 'the track can be coloured by any of these three variables';
echo 'i.e. -c CELL would colour all tracks from a cell line the same colour.';
echo '-c CONDITION would colour all conditions (e.g. H3K4me3..) the same';
fi
shift $(($OPTIND -1))