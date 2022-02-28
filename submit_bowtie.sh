#!/bin/sh

#submit script for run_Bowtie_on_VACC.sh
#OPTIONS
#ADD OPTARGS

while getopts "i:g:hw::nt5::3::f" opt; do
    case $opt in
        i) FASTQFOLDER=$OPTARG;;
        g) GENOME=$OPTARG;;
        w) WD=$OPTARG;;
        n) FNC="Y";;
        t) TRIM="Y";;
        5) TRIM5=$OPTARG;;
        3) TRIM3=$OPTARG;;
        f) FILTER="Y";;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
    esac
done
shift $(($OPTIND -1))

if [[ -z $GENOME ]]; then echo "GENOME (-g) is not set"; HELP=1; else echo "The following genome will be used: $GENOME"; fi
if [[ -z $FASTQFOLDER ]]; then echo "INPUT DIRECTORY (-i) is not set"; HELP=1; else echo "The path to fastq files is: $FASTQFOLDER"; fi
if [[ -z $WD ]]; then echo ""; else echo "The working directory is: $WD"; OTHEROPTS=",WD=$WD"; fi
if [[ -z $FNC ]]; then echo "File naming convenion not set"; else echo "File naming convenion is set"; OTHEROPTS="$OTHEROPTS,FNC=$FNC"; fi
if [[ -z $TRIM ]]; then echo ""; else echo "Reads will be trimmed before alignment"; OTHEROPTS="$OTHEROPTS,FILTER=$TRIM"; fi
if [[ -z $TRIM5 ]]; then echo ""; else echo "5\' ends of reads will be trimmed by $TRIM5 bases"; OTHEROPTS="$OTHEROPTS,TRIM5=$TRIM5"; fi
if [[ -z $TRIM3 ]]; then echo ""; else echo "3\' ends of reads will be trimmed by $TRIM3 bases"; OTHEROPTS="$OTHEROPTS,TRIM3=$TRIM3"; fi
if [[ -z $FILTER ]]; then echo ""; else echo "READS will be filtered by TRIMMOMATIC"; OTHEROPTS="$OTHEROPTS,FILTER=$FILTER"; fi
if [[ $HELP == 1 ]]; then
echo 'SUBMITS A BATCH OF BOWTIE JOBS';
echo 'REQUIRES (-g) A genome: either "mm10" or "hg19"';
echo 'REQUIRES (-i) A path to the folder containing fastq files';
echo 'OPTION: (-w) A working directory to output results. If this is not set th bowtie script will create one';
echo 'OPTION: (-n) If flagged it will assume that files are in CELL_CONDITION_REPLICATE format for file naming';
echo 'OPTION: (-t) If flagged this will trim the fastq reads using Trimmomatic before aligning';
echo 'OPTION: (-5) If flagged this will trim the 5 prime ends of fastq in bowtie by x base pairs. Requires a value';
echo 'OPTION: (-3) If flagged this will trim the 3 prime ends of fastq in bowtie by x base pairs. Requires a value';
echo 'OPTION; (-h) This information';
exit 1;
else echo ""; fi

# CHECK THE FILES HERE
FASTQS="$FASTQFOLDER/*.fastq*"


OPTIONS=",GENOME=$GENOME""$OTHEROPTS"
echo "$OPTIONS"
# COUNTER
TEMPFILE=~/.tmp; echo 0 > $TEMPFILE;

for f in $FASTQS; do
#LOOP COUNTER
COUNTER=$[$(cat $TEMPFILE) + 1]; echo $COUNTER > $TEMPFILE
NAME="BOWTIE""$COUNTER"

qsub -v INPUT=$f"$OPTIONS" -N $NAME /users/j/a/jargordo/bin/scripts/run_bowtie_on_VACC.sh
echo " -v parameters: qsub -v INPUT=$f""$OPTIONS -N $NAME"
done
unlink $TEMPFILE
