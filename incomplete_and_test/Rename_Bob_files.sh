#!/bin/sh

#converts files from Bioniformatics/Bob Devins to something more manageable
#REQUIRES FILE PATHS (FOLDER, MVFOLDER)
# The FOLDER should look something like this: /users/r/d/rdevins/scratch/MPS/140225_SNL128_0064_AC30CGACXX/samples_out


while getopts "f:d:h" opt; do
    case $opt in
        f) FOLDER=$OPTARG;;
        d) MVFOLDER=$OPTARG;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))

if [ -z $FOLDER ]; then "THE FOLDER TO COPY FROM WAS NOT SET OR DOES NOT EXIST" exit 1; else echo "FOLDER WAS SET"; fi
if [ -z $MVFOLDER ]; then MVFOLDER="/users/j/a/jargordo/scratch" echo "MVFOLDER WAS NOT SET: NOW SET TO DEFAULT"; else echo "MVFOLDER WAS SET"; fi
echo -e "----------------------------------------------------------------------\n\n"

cd $FOLDER
#Inside the FOLDER the are several files and folders you need:
#There should be a demultiplex file named "demux_*.csv". The * should be the lane number that your sample ran on.
#This is the demultiplex file for barcode assignment. You should double check this file to make sure that the barcodes are correct for the sample.
#IN FACT it is so important we will just cat it here:
DEMUXFILE="../"$FOLDER/"*.csv"
echo "SAMPLES WERE DEMULTIPLEXED USING THIS SCHEME:"
cat $DEMUXFILE

#You don't really need anything else from this parent FOLDER unless things go really wrong.
#Move to LANE specific folder "Unaligned_Lane*"
LANEFOLDER="sample_out"

#The number of reads after Basecall/demultiplexing are located in a html file in the Basecall_Stats_* folder
#copy that over
cp ./$LANEFOLDER/Basecall_Stats_*/Demultiplex_Stats.htm "$MVFOLDER/Demultiplex_Stats_""$LANEFOLDER"".htm"

#The important files (fastq) are in the Project folder that looks something like this "Project_jstein7_bccl_design_20131022_Lane1"
PROJECT=$FOLDER/$LANEFOLDER/Project*/*
#should be multiple folders matching the number of samples/barcodes
#These folders are usually named "Sample_Name_Name"

for f in $PROJECT; do
PROJECTFOLDER=$(echo $f | awk -F '/' '{print $NF}')
#Remove the "Sample" and convert the name to uppercase (to make the greps easier)
LSAMPLE=$(echo $PROJECTFOLDER | awk -F 'Sample_' '{$1=""; print toupper ($0)}')

#Sample should be identified by CELL_CONDITION_REPLICATE_...
#However that doesn't always happen so... we will search for a HISTONE (H3K or H4K), a known CELL and if DAY appears somewhere
#Also look for a R*, "REP" or the word "REPLICATE"
HISTONE=$(echo $LSAMPLE | grep -Eo 'H3K|H4K[^_]*')
CELL=$(echo $LSAMPLE | grep -Eo 'MDA*|MCF*|MCF_7|MCF_10*|BMSC|MC3*[^_]*')
DAY=$(echo $LSAMPLE | grep -Eo 'DAY[0-99]|DAY_[0-99][^_]*'); DAY=${DAY//_/}; DAY=$(echo $DAY | sed  's/DAY\([0-9]\)/DAY0\1/');
REPLICATE=$(echo $LSAMPLE | grep -Eo 'R[1-9]|REP[1-9]|REP_[1-9]|REPLICATE[1-9]|REPLICATE_[1-9][^_]*');REPLICATE=${REPLICATE//_/};

#Clean up the naming
if [ $CELL = "MDA" ]; then CELL="MDA231";
    elif [ $CELL1 = "MCF_7" ]; then CELL="MCF7";
    elif [ $CELL1 = "MCF_10A" ]; then CELL="MCF10A";
else echo ""; fi

#Check if any of these variables were assigned
FNAME=""
TEMP=~/.tmp; echo 0 > $TEMP;
if [ -z $CELL ]; then echo -e "NO CELL LINE WAS IDENTIFIED OR SET IN THE QSUB\nUSING THE PLACE HOLDER \"SOMECELL\"\n"; FNAME="SOMECELL"; else FNAME="$CELL"; fi
if [ -z $HISTONE ]; then COUNT=$[$(cat $TEMP) + 1]; echo $COUNT > $TEMP; else FNAME=$FNAME"_"$HISTONE; fi
if [ -z $DAY ]; then COUNT=$[$(cat $TEMP) + 1]; echo $COUNT > $TEMP; else FNAME=$FNAME"_"$DAY; fi
if [ -z $REPLICATE ]; then COUNT=$[$(cat $TEMP) + 1]; echo $COUNT > $TEMP; else FNAME=$FNAME"_"$REPLICATE; fi
COUNT=$(cat $TEMP); if [ $COUNT < 2 ]; then FNAME=$FNAME"_"$LSAMPLE; else echo ""; fi
unlink $TEMP
#Find barcodes from the DEMUXFILE and append them to the file name to make it a little easier to troubleshoot
BARCODE=$(grep -o -P "(?<=$LESSSAMPLE,human|mouse,).*(?=,RNaseq|ChIPseq)" $DEMUXFILE)


FNAME1=$CELL"_"$LESSSAMPLE"_"$BARCODE"_"$OTHERDETS
FASTQC="$MVFOLDER"/"$FNAME1""$FNAME2"
FASTQ=$MVFOLDER/$FNAME1".fastq.gz"

#check for Concat21
cp $f/Concat11_fastqc.zip $FASTQC
cp $f/Concat11.fastq.gz $FASTQ
done

echo "THE FILES WRE COPIED FROM: $FOLDER"
echo "THE FILES WERE COPIED TO: $MVFOLDER"

