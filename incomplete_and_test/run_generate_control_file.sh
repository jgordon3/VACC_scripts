#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash


# Accept file prefix from run_MACS2_on_SGE.sh (and other scripts to generates a UCSC "control" file for bigwigs and bigbeds
# qsub -V TREATABRVNAME= ~/scripts/run_generate_control_file.sh

#######################################################
####################### DEFAULTS GO HERE ##############
#######################################################
URL=https://galaxy.med.uvm.edu/static/UCSCtracks/hg38/ChIP_seq/
URL_FOLDER=/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/hg38/ChIP_seq/
VIEW_LIMITS=0.001:3
SMOOTH_WIN=2
MAX_HEIGHT=128:24:8
COLOR="166,206,227"
VIS="full"

usage ()
{
DIV1=`eval printf '=%.0s' {1..100}`; DIV2=`eval printf '=%.0s' {1..25}`
printf "%s\n" "" "$DIV1"
printf "%s\t" "$DIV2" "run_generate_control_file.sh" "" "" "" "$DIV2"
printf "%s\n" "" "$DIV1" ""
printf "%s\n" "Accept file prefix from run_MACS2_on_SGE.sh (and other scripts to generates a UCSC control file for bigwigs and bigbeds"
printf "%s\n" "qsub -V TREATABRVNAME= ~/scripts/run_generate_control_file.sh"
printf "%s\n" "" "$DIV1"
printf "%s\n"
exit 1;
}

# OPTIONS
while getopts "T:i:c:g:p:m:u:s:a:h:b:w:q:v:" opt; do
    case $opt in
        "-T"|"--treatment"  ) TREAT=$OPTARG;;          #treatment name is needed"
        "-i"| "--input" ) CONTROL=$OPTARG;;
        "-c"|"--compress"   ) COMPRESS=$OPTARG;;
        "-g"| "--genome"    ) GENOME=$OPTARG;;
        "-p"|"--pvalue_cutoff"  ) PVALUE=$OPTARG;;
        "-m"|"--mfold_cutoff"   ) MFOLD=$OPTARG;;
        "-u"|"--user_name"  ) NAME=$OPTARG;;
        "-s"|"--save_bedgraph"  ) SAVEBED=$OPTARG;; #save bedgraph output from MACS or peakcaller
        "-a"|"--add_peakcall"   ) ADD_PEAKCALL=$OPTARG;;
        "-h"|"--help"   ) usage;;  #runs help file
        "-w"|"--URL" ) URL=$OPTARG;; #URL address for UCSCGB access
        "-q"|"--URL_FOLDER" ) URL_FOLDER=$OPTARG;; #web accessable folder to transfer tracks
        "-v"|"-veiw_limits" ) VIEW_LIMITS=$OPTARG;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))


#if [[ ! -e $TREAT ]]; then echo "Need valid path to BOWTIE BAM/SAM file"; usage; fi
#if [[ ! -e $CONTROL ]]; then echo "No Input set"; fi


# give flag option to just assign all bigWigs in the folder (*.bw)
# if [[ ! -z $ALL ]]; then

# check/conditional these later in case script is launched independently
#Check if files exist; if not check for tar archive

TREATWIG=$TREATABRVNAME"_sorted_CPM_norm.bw"
CTRLWIG=$CONTROLABRVNAME"_sorted_CPM_norm.bw"
TREATFEWIG=$TREATABRVNAME"_sorted_FE.bw"
TREATlogLRWIG=$TREATABRVNAME"_sorted_logLR.bw"
TREATBIGBED=$TREATABRVNAME"_sorted_narrowPeaks.bb"



#check to see if this exists and delete it
FILE=$TREATABRVNAME"_UCSC_control_file.txt"
if [[ -e $FILE ]]; then rm $FILE; fi

for f in $TREATWIG $CTRLWIG $TREATFEWIG $TREATlogLRWIG; do
    printf "track type=bigWig name=$f description="'"'"$f"'"' >> $FILE
    printf " bigDataUrl=$URL/$f autoScale=off viewLimits=$VIEW_LIMITS smoothingWindow=$SMOOTH_WIN windowingFunction=maximum maxHeightPixels=$MAX_HEIGHT color=$COLOR visibility=$VIS \n" >> $FILE
done

for i in $TREATBIGBED; do
printf "track type=bigBed name=$i description="'"'"$i"'"' >> $FILE
printf " bigDataUrl=$URL/$i visibility=$VIS \n" >> $FILE
done

#move everything to URL folder

####### HELP #########




