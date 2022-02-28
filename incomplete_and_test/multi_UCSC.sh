#!/bin/bash

# DEFAULTS GO HERE

URL="https://galaxy.med.uvm.edu/static/UCSCtracks/Galaxy_generated/"; HOST="/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/Galaxy_generated";
GENOME="db=hg38";CHR_POSITION=""; VIEW_LIMITS=0:10; SMOOTH_WIN=2; MAX_HEIGHT=128:24:8
VIS="full"

####### HELP #########

usage ()
{
printf "%s\n" "UCSC OPERATIONS"
printf "%s\n" "This script is called from the macs2.sh script to generate UCSC visualizations"
printf "%s\n" "REQUIRES (NON OPTIONAL): A folder containing bigwigs and bigbed files assigned with the -F flag"
printf "%s\n" "Looks for all bigwig (*.bw) and bigbeds (*.bb) and wil move and generate a UCSC control file"
printf "%s\n" "OPTIONS: -U assign a different URL to host UCSC tracks (DEFAULT: $URL)"
printf "%s\n" "OPTIONS: -H assign a different folder to host UCSC tracks (DEFAULT: $HOST)"
printf "%s\n" "OPTIONS: -G assign a different genome for UCSC tracks (DEFAULT: $GENOME)"
printf "%s\n" "OPTIONS: -p assign a different starting postion for UCSC tracks (DEFAULT: $CHR_POSITION)"
exit 1;
}


# OPTIONS
while getopts "T:U:H:p:h" opt; do
    case $opt in
        F) FOLDER=$OPTARG;;
        U) URL=$OPTARG;;
        H) HOST=$OPTARG;;
        G) GENOME=$OPTARG;;
        p) CHR_POSITION=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

#assign a random color if variable is not passed
#RED/ORANGE
C1="169,58,3";C2="190,137,10";C3="114,37,15";C4="179,46,43";C5="236,79,21";C6="200,8,21";
#BLUE/GREEN
C7="18,73,93"; C8="71,99,136";C9="120,143,83";C10="60,62,26";C11="17,61,28";C12="25,28,47";C13="13,37,76";
N="C"$(( ( RANDOM % 12 )  + 1 ))
RANDOMCOLOR=${!N}
if [ -z "$COLOR" ]; then COLOR="$RANDOMCOLOR"; fi


# Generate UCSC file
UCSC_FILE=`basename $FOLDER`
UCSC_FILE="${UCSC_FILE}_ucsc_file.txt"

# bigwigs

for f in $FOLDER/*.bw $FOLDER/*.bigwig; do
    TRACK_NAME=`basename $f`
    TRACK_NAME=${TRACK_NAME%.*}
	printf "track type=bigWig"
    printf "%s " "name=\"$TRACK_NAME\"" "description=\"galaxy generated bigwig for $TRACK_NAME\"" "bigDataUrl=$URL$f" "autoScale=off" "viewLimits=$VIEW_LIMITS" >> $UCSC_FILE
	printf "%s " "smoothingWindow=$SMOOTH_WIN" "windowingFunction=maximum" "maxHeightPixels=$MAX_HEIGHT" "color=$COLOR" "visibility=$VIS" "\n" >> $UCSC_FILE

done

for f in $FOLDER/*.bb $FOLDER/*.bigbed; do
    TRACK_NAME=`basename $f`
    TRACK_NAME=${TRACK_NAME%.*}
	printf "track type=bigBed" >> $UCSC_FILE
	printf "%s " "name=\"$TRACK_NAME\"" "description=\"galaxy generated bed for $track_name\"" "bigDataUrl=$URL$track" >> $UCSC_FILE

	printf " visibility=$VIS" >> $FILE
	printf "\n" >> $FILE
done

# html file 

printf "https://genome.ucsc.edu/cgi-bin/hgTracks?$GENOME&position=chr7:27077339-27212428&hgct_customText=$URL$FILE" >> $HTML_FILE

#copy the control file to folder
cp $FILE $URL_FOLDER



#move everything to URL folder
