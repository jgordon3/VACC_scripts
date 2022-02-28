#!/bin/bash

# run_demux.sh
# bcl2fastq conversion for files from AGTC
# no sanity checks

source /slipstream/home/jonathan/scripts/qstat_query.sh

while getopts "s:d:l:h" opt; do
    case $opt in
        s) SAMPLE=$OPTARG;;#sample sheet
        l) LANE=$OPTARG;;
        d) DATAFOLDER=$OPTARG;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
if [ -z $SAMPLE ]; then echo "Need valid path to sample sheet"; exit 1; fi
if [ -z $DATAFOLDER ]; then echo 'Need a full path to the data folder'; exit 1; else echo "Data folder: $DATAFOLDER"; cd $DATAFOLDER; fi
shift $(($OPTIND -1))


# checks for Data folder
if [ -d ./Data/Intensities ]; then echo "Found data ==> starting configuration";
    else echo "This folder does not appear to have Data to convert. Try relaunching script from a folder containing data or using the -d flag."; exit 1;
fi

# count number of lanes in DATAFOLDER
NUM_FOLDERS=$(find $DATAFOLDER/Data/Intensities/L00* -maxdepth 0 -type d | wc -l)
echo $NUM_FOLDERS


BASELANE=`basename $LANE`
OUTPUTDIR="$BASELANE""_demultiplex"

# loop through data folder
#i=1
#while [ $i -le $NUM_FOLDERS ]; do
#    HOLDER=L00$i
#    HOLDER2=$HOLDER"_demultiplex"
#    echo $HOLDER2
#    i=$((i+1))
#done

#  =========== SAMPLESHEET OPERATIONS =============
# UNDER CONSTRUCTION
# want to run lane by lane (seperate qsubs) without opening the samplesheet
# convert the AGTC sample sheet to a usable version

# convert to true csv from mac style lines
perl -w015l12pi -e1 $SAMPLE

# split samplesheet
#SAMPLEBASE=${SAMPLE%.csv}
#OLDIFS=$IFS
#IFS=","
#REMOVE FIRST LINE
#while read FCID Lane SampleID SampleRef Index Description Control Recipe Operator SampleProject || [ "$Lane" ]; do
#        SS="$SAMPLEBASE""_""L00$Lane"
#        echo "$FCID $Lane $SampleID $Index"  >> $SS
#done < $SAMPLE
#IFS=$OLDIFS
# ===================================================================

# configure the bcl2fastq
/usr/local/bin/configureBclToFastq.pl --input-dir ./Data/Intensities/BaseCalls --output-dir $OUTPUTDIR --sample-sheet $SAMPLE --no-eamss --fastq-cluster-count 0

sleep 5
cd $OUTPUTDIR

#copy sample sheet into project folder
cp $SAMPLE $OUTPUTDIR
# make wrapper and qsub
cp /slipstream/home/jonathan/scripts/make_wrapper.sh ./
# write -o line (sge_job folder) to current directory

JOB1=$(qsub make_wrapper.sh)
TRIM_JOBID=$(awk -v RS=[0-9]+ '{print RT+0;exit}' <<< "$JOB1") #returns JOBID
TRIM_JOBS="$TRIM_JOBS $TRIM_JOBID"

qstat_query($TRIM_JOBS)

BASECALL_STATS="$OUTPUTDIR""/Basecall_Stats_*/Demultiplex_Stats.htm"
BASECALL_OUT="$OUTPUTDIR""_stats.htm"
cp ./$BASECALL_STATS ./$BASECALL_OUT


# ======= TO DO ===================
# Copy fastq.gzs, tar and move
# Archive folder (tar.gz)
# write to galaxy






