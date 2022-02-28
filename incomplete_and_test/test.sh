#!/bin/sh

#!/bin/sh

# GENERATES CONTROL FILE FOR UCSC GENOME BROWSER FROM ALL* FILES IN A DIRECTORY
# REQUIRES: Control file name (-n). Program will remove file if it exists.
# REQUIRES:
# OPTIONS: validate files before writing:

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
shift $(($OPTIND -1))

if [ -z $NAME ]; then echo "Need a name for the control file (-n)"; exit 1; else echo "$NAME"; fi
if [[ ! -z $VAL ]]; then echo "Files will be validated before writing"; fi
if [ -z $HELP ]; then echo ""; else
echo 'GENERATES A UCSC GENOME BROWSER CONTROL FILE FROM EVERY FILE IN A DIRECTORY'
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
echo '-f [FILEPATH] to directory with files. Default is current directory'; fi
if [ -z $FILEPATH ]; then FILEPATH=$(pwd); echo "FILEPATH is set to the current directory"; else echo "FILEPATH was set to $FILEPATH"; fi
if [ -z $URL ]; then URL2=$(echo "$(pwd)" | awk -F 'UCSC_data_feed' '{print $2"/" }');  URL="http://www.med.uvm.edu/steinlab/""$URL2"; echo "URL is set to the current directory"; else echo "URL was set to $URL"; fi
if [ -f "$NAME" ]; then echo "$NAME exists: removing now"; rm $NAME; else echo "writing to $NAME"; fi

FILES=$FILEPATH"/*"
echo "$FILES"


#COLOR BREWER 12 COLOR PALETTE
C1="166,206,227";C2="31,120,180";C3="178,223,138";C4="51,160,44";C5="251,154,153";C6="227,26,28";C7="253,191,111";
C8="255,127,0";C9="202,178,214";C12="106,61,154";C11="255,255,153";C12="177,89,40";

APREV=0
CELLPREV=0
i=1
CPREV=$(echo "C" $i)
echo "$CSEP"

for f in $FILES; do
# let me know it is working...
echo "Processing $f file..."

#stupid (but quick) way of breaking down the file name
FILE=$(echo "$f" | awk -F '/' '{print ( $(NF) ) }')
CONDITION=$(echo "$FILE" | awk -F '_' '{print $2}'| awk '{print toupper($0)}')
CELL=$(echo "$FILE" | awk -F '_' '{print $1}' | awk '{print toupper($0)}')
REPLICATE=$(echo "$FILE" | awk -F '_' '{print $3}')
LONGDESC=$(echo "$FILE" | awk -F '.' '{print $1}')
FILETYPE=$(echo "$FILE" | awk -F '.' '{print $2}')

# SET COLOURS

if [ $CSEP = "CELL" ]; then
    if [ "$CELL" = "MCF7" ]; then COLOUR=$C1; echo $COLOUR;
        elif [ "$CELL" = "MDA231" ]; then COLOUR=$C3; echo $COLOUR;
        elif [ "$CELL" = "MCF10A" ]; then COLOUR=$C5; echo $COLOUR;
        elif [ "$CELL" = "BMSC" ]; then COLOUR=$C7; echo $COLOUR;
        else echo "$CELL not in the normal list";
            if [ "$CELL" = "$CELLPREV"]; then COLOUR=$CPREV;
            else i=$((i+1)); CPREV=$(echo "C" $i); COLOUR=$CPREV;
            fi
    fi
elif [ $CSEP = "CONDITION" ]; then echo "set to condition";
else echo "CSEP is something else";
fi
done