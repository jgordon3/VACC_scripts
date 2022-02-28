#!/bin/sh

#Remote access
#ssh $VACC; cd scratch;
#OPTIONS

while getopts "i:s:w:hr:" opt; do
    case $opt in
        i) FASTQFOLDER=$OPTARG;;
        s) SPECIES=$OPTARG;;
        w) WD=$OPTARG;;
        r) REFERENCE=$OPTARG;;
        n) FNC="Y";;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))
if [[ -z $SPECIES ]]; then echo "SPECIES (-s) is not set"; HELP=1; else echo "The following genome will be used: $SPECIES"; fi
if [[ -z $FASTQFOLDER ]]; then echo "INPUT DIRECTORY (-i) is not set"; HELP=1; else echo "The path to fastq files is: $FASTQFOLDER"; fi
if [[ -z $WD ]]; then echo 'The working directory (-w) is not set. Setting a default "TOPHAT_DATE.." directory';
    else echo "The working directory is: $WD"; fi
if [[ -z $REFERENCE ]]; then echo "REFERENCE (-r) is set to UCSC by default"; else echo "The REFERENCE was set to: $REFERNCE by -r"; fi
if [[ $HELP == 1 ]]; then
echo 'SUBMITS A BATCH OF TOPHAT JOBS through the script run_TOPHAT_on_VACC.sh ';
echo 'REQUIRES (-s) A species: either "mm10" or "hg19"';
echo 'REQUIRES (-i) A path to the folder containing fastq files';
echo 'OPTION (-r) A reference to align to "NCBI" or "UCSC". Default is UCSC';
echo 'OPTION: (-w) A working directory to output results. If this is not set the script will create one';
echo 'OPTION:';
echo 'OPTION; (-h) This information';
exit 1;
else echo ""; fi

FASTQS=$FASTQFOLDER"/*"

# CHECK FILE NAMES
if [[ -z $FNC ]]; then echo "File naming convenion not set";


else echo "File naming convenion is set"; fi


#while true; do
#read -p "Do you wish to install this program?" yn
#case $yn in
#[Yy]* ) make install; break;;
#[Nn]* ) exit;;
#* ) echo "Please answer yes or no.";;
#esac
#done


SCRIPT=/users/j/a/jargordo/bin/scripts/run_TOPHAT_on_VACC.sh
# ARE FILES IN THE FOLLOWING FORMAT: "CELL_CONDITION_RELICATE_..." (Y or N)
FNC="Y"

#1) USE A WHOLE FOLDER


# COUNTER
TEMPFILE=~/.tmp; echo 0 > $TEMPFILE;

for f in $FASTQS; do
#LOOP COUNTER
COUNTER=$[$(cat $TEMPFILE) + 1]; echo $COUNTER > $TEMPFILE
NAME="TOPHAT""$COUNTER"

qsub -v PATHTOFASTQ=$f,SPECIES=$SPECIES,REFERENCE=$REFERENCE,FNC=$FNC,WD=$WD,LIBTYPE=$LIBTYPE -N $NAME $SCRIPT
done
unlink $TEMPFILE
