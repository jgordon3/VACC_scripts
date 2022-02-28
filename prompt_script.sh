#!/bin/sh
#  prompt_script.sh
# This will prompt the user for a series of options to submit a TOPHAT job(s)
# This script writes the variables straight to a copied run_TOPHAT_on_VACC.sh script
PARENTSCRIPT=/users/j/a/jargordo/bin/scripts/run_TOPHAT_on_VACC.sh
CUFFLINKSCRIPT=/users/j/a/jargordo/bin/scripts/run_CUFFLINKS_on_VACC.sh

printf "\n"
printf "This script will submit the TOPHAT script to analyze RNA-seq fastq data. \n"
printf "There will be a series of prompts. Please provide answers. \n"
printf "The entire RNA-SEQ pipeline can be launched from this script (if selected). \n"
printf "It will also write scripts for individual parts (TOPHAT, CUFFLINKS, CUFFDIFF) \n"
printf "Please use absolute paths when providing locations. \n"
printf "\n"

#ASSIGN PATHTOFASTQS  (CHECK)
read -e -p "Path to folder containing fastq file(s): " PATHTOFASTQS
#PATHTOFASTQS=$(echo $PATHTOFASTQS | sed 's/ /\\ /')
PATHTOFASTQS=${PATHTOFASTQS%/}
eval PATHTOFASTQS=$PATHTOFASTQS
if [ -d $PATHTOFASTQS ];
    then cd "$PATHTOFASTQS"; FILES=$PATHTOFASTQS"/*"; printf "\nThis appears to be a directory...checking for fastq files: \n"
    elif [ -f "$PATHTOFASTQS" ]; then printf "\nEntered one file: \n"; FILES=$PATHTOFASTQS;
fi

#COUNT FILES
GZCOUNT=0;FASTQCOUNT=0;OTHERCOUNT=0;
for f in $FILES; do
EXT=$(echo $f | awk -F "." '{print toupper($NF)}')
    if [[ "$EXT" = "GZ" ]]; then EXT2=$(echo $f | awk -F "." '{print toupper($(NF-1))}');
        if [[ "$EXT2" = "FASTQ" ]]; then GZCOUNT=$((GZCOUNT+1)); FASTQFILES="$FASTQFILES $f"; fi
    elif [[ "$EXT" = "FASTQ" ]]; then FASTQCOUNT=$((FASTQCOUNT+1)); FASTQFILES="$FASTQFILES  $f";
    else OTHERCOUNT=$((OTHERCOUNT+1));
    fi
done
SUBMITCOUNT=$((GZCOUNT+FASTQCOUNT));
#if less than 0 exit
printf "There is/are $FASTQCOUNT fastq file(s) to process: \n"
if [[ $GZCOUNT -ge 1 ]]; then printf "There are $GZCOUNT file(s) that appear to be gzipped fastq file(s). \n"; fi
if [[ $OTHERCOUNT -ge 1 ]]; then printf "There are $OTHERCOUNT file(s) are neither fastq or gziped fastq (fastq.gz). \n"; fi
if [[ $SUBMITCOUNT -lt 1 ]]; then printf "It does not look like there are any fastqs to process. Please check extensions (.fastq, .fastq.gz) or choose a new folder. \n"; exit; fi
printf "\n"

#PAIRED END
while true; do
    read -p "Is this data paired-end: (y or n)? " yn
        case $yn in
        [Yy]* ) PAIRCHECK=1; break;;
        [Nn]* ) PAIRCHECK=0; break;;
        * ) echo "\n"; echo "Please answer yes or no.";
    esac
done
#FIND PAIRS
if [ $PAIRCHECK == 1 ]; then
printf "Attempting to find pairs: \n"
printf "\n"
printf "Paired files are: \n"
echo $FASTQFILES | xargs ls -l | sort -k 5,6 |
    while read PERMS LINKS USER GROUP SIZE M D Y FILE; do
        if [[ "$SIZE" -eq "$LASTSIZE" ]]; then
            F1=$(echo $LASTFILE | awk -F "/" '{print $NF }'); F2=$(echo $FILE | awk -F "/" '{print $NF }')
            echo "$F1 and $F2"; echo "$LASTFILE $FILE" >> pairs.txt
        else
            LASTSIZE="$SIZE"; LASTFILE="$FILE"
        fi
    done
fi
printf "\n"
#is this correct?
#Do uo want to submit all these files?

#ASSIGN WD
read -e -p "Working directory for output: " PATHTOWD
PATHTOWD=${PATHTOWD%/}
eval PATHTOWD=$PATHTOWD
if [ -d $PATHTOWD ]; then cd $PATHTOWD; else echo "Directory does not exist... creating..."; mkdir $PATHTOWD; cd $PATHTOWD; fi
printf "\n"

#ASSIGN SPECIES
while true; do
    read -p "Which genome do you want to align to? Options are: mm10 or hg19 (m or h):" mh
        case $mh in
        [Mm]* ) SPECIES="mm10"; break;;
        [Hh]* ) SPECIES="hg19"; break;;
        * ) echo "\n"; echo "Please select m (mm10) or h (hg19).";
    esac
done
printf "\n";

#ASSIGN REFERENCE
while true; do
    read -p "Which reference do you want to use for the GTF file? Options are: UCSC (u), NCBI (n), ENSEMBL (e) or GENCODE (g):" uneg
        case $uneg in
        [Uu]* ) REFERENCE="UCSC"; break;;
        [Nn]* ) REFERENCE="NCBI"; break;;
        [Ee]* ) REFERENCE="ENSEMBL"; break;;
        [Gg]* ) REFERENCE="GENCODE"; break;;
        * ) echo "\n"; echo "Please select u (UCSC), n (NCBI), e (ENSEMBL) or g (GENCODE).";
    esac
done
printf "\n";
printf "You have enough information to submit the job with the default parameters: \n";
printf "$SUBMITCOUNT file(s) to submit, which is(are): \n";
for f in $FASTQFILES; do printf "$f \n"; done;
printf "\n";
printf "Reads will be aligned to $SPECIES genome using $REFERENCE gene annotation \n";
printf "\n";

#CONTINUE CHECK
while true; do
    read -p "Do you want to proceed to advanced options (p), submit (s), write scripts and exit (w), or abort (q)? " pswq
        case $pswq in
        [Pp]* ) ADVANCED=1; break;;
        [Ss]* ) ADVANCED=0; break;;
        [Ww]* ) ADVANCED=2; break;;
        [Qq]* ) printf "\n"; echo "Submission aborted"; exit; break;;
        * ) echo "\n"; echo "Please answer proceed (p), submit (s) or quit (q).";
    esac
done
printf "\n";

#START CREATING SCRIPTS HERE: write variables to script for each file
for f in $FASTQFILES; do

    FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
    SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
    JOBNAME="#PBS -N ""$FNAME""_TOPHAT"; awk -v a="$JOBNAME" '{ if(NR==7) { print a} else {print $0} } ' $PARENTSCRIPT > $SHSCRIPT;
    awk -v b="PATHTOFASTQ=$f;" '{ if(NR==14) { print b} else {print $0} }' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    awk -v c="WD=$PATHTOWD;" '{ if(NR==15) { print c} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    awk -v d="SPECIES=$SPECIES;" '{ if(NR==16) { print d} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    awk -v e="REFERENCE=$REFERENCE;" '{ if(NR==17) { print e} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    chmod a+x $SHSCRIPT
done

#DEFAULT PARAMETERS EXIT
if [ $ADVANCED == 0  ]; then
    echo `qsub $SHSCRIPT`; exit;
elif [ $ADVANCED == 2  ]; then
    printf "TOPHAT scripts were written to the working directory \n"; exit;
else continue;
fi

#TRIM READS
while true; do
    read -p "Do you want trim reads for quality and/or length using TRIMMOMATIC before alignment: (y or n)? " yn
        case $yn in
        [Yy]* ) TRIMCHECK=1; break;;
        [Nn]* ) TRIMCHECK=0; break;;
        * ) echo "\n"; echo "Please answer yes or no.";
    esac
done
if [ $TRIMCHECK == 1 ]; then
read -e -p "Please enter a number for quality cutoff (Default is 20) : " SLIDINGWINDOW
read -e -p "Please enter a number for trimming 5' Ends (Default is 3) : " LEADING
read -e -p "Please enter a number for trimming 3' Ends (Default is 3) : " TRAILING
SLIDINGWINDOW="4:$SLIDINGWINDOW"
#add null check here
    for f in $FASTQFILES; do
        FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
        SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
        awk -v s="SLIDINGWINDOW=$SLIDINGWINDOW;" -v l="LEADING=$LEADING;" -v t="TRAILING=$TRAILING;" -v x="TRIMCHECK=$TRIMCHECK;" '{ if(NR==18) { print x,s,l,t} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"
        mv "$SHSCRIPT.temp" $SHSCRIPT
    done
else continue;
fi

#STRAND TYPE
while true; do
    read -p "Do you want set the strand type (--library-type) (default is unstranded): (y or n)? " yn
        case $yn in
        [Yy]* ) STRCHECK=1; break;;
        [Nn]* ) STRCHECK=0; break;;
    * ) echo "\n"; echo "Please answer yes or no.";
    esac
done
if [ $STRCHECK == 1 ]; then
    while true; do
        read -p "Set strand type (--library-type): options are: fr-unstranded (u), fr-firststrand (f) or fr-secondstrand (s) " ufs
            case $ufs in
            [Uu]* ) LIBTYPE="fr-unstranded"; break;;
            [Ff]* ) LIBTYPE="fr-firststrand"; break;;
            [Ss]* ) LIBTYPE="fr-secondstrand"; break;;
            * ) echo "\n"; echo "Please select u (fr-unstranded), f (fr-firststrand) or s (fr-secondstrand). \n"; echo "If you are unsure consult the TOPHAT manual. The default is fr-unstranded (u) \n"
        esac
    done
for f in $FASTQFILES; do
    FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
    SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
    awk -v f="LIBTYPE=$LIBTYPE;" '{ if(NR==19) { print f} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
done
else continue;
fi

#ANCHOR LENGTH
while true; do
    read -p "Do you want set the anchor length (-a) (default is 8): (y or n)? " yn
        case $yn in
        [Yy]* ) ALENCHECK=1; break;;
        [Nn]* ) ALENCHECK=0; break;;
        * ) echo "\n"; echo "Please answer yes or no. ";
    esac
done
if [ $ALENCHECK == 1 ]; then
    read -e -p "Please enter a number for anchor length: " ALENGTH
    eval ALENGTH=$ALENGTH
    for f in $FASTQFILES; do

    FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
    SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
    awk -v g="ALENGTH=$ALENGTH;" '{ if(NR==20) {print g} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    done
else continue;
fi

#JUNCTIONS
while true; do
    read -p "Do you want to alter junction handling: (y or n)? " yn
        case $yn in
        [Yy]* ) JUNCHECK=1; break;;
        [Nn]* ) JUNCHECK=0; break;;
        * ) echo "\n"; echo "Please answer yes or no.";
    esac
done
if [ $JUNCHECK == 1 ]; then
    while true; do
        read -p "Junction options: --no-novel-juncs (n), --no-gtf-juncs (g) or use a supplied junction file with no novel junctions (j) " ngj
            case $ngj in
            [Nn]* ) JUNCTIONS="--no-novel-juncs"; break;;
            [Gg]* ) JUNCTIONS="--no-gtf-juncs"; break;;
            [Jj]* ) JUNCTIONS="--no-novel-juncs";
                    read -e -p "Please enter a filepath for the juncs file: " JUNCTIONF
#check this later
                    break;;
            * ) echo "\n"; echo "Please select n (no novel junctions), g (no gtf junctions) or j (supply a junction file). \n"; echo "If you are unsure consult the TOPHAT manual. \n"
        esac
    done
    for f in $FASTQFILES; do
        FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
        SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
        awk -v h="JUNCTIONS=$JUNCTIONS;" '{ if(NR==21) {print h} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
        if [ ! -z $JUNCTIONF ]; then
        awk -v h="JUNCTIONF=$JUNCTIONF;" '{ if(NR==22) {print h} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
        fi
    done
else continue;
fi

#BIGWIGS
while true; do
    read -p "Do you want to create bigwig files to display on UCSC: (y or n)? " yn
        case $yn in
        [Yy]* ) BAMCHECK=1; break;;
        [Nn]* ) BAMCHECK=0; break;;
        * ) echo "\n"; echo "Please answer yes or no.";
    esac
done
if [ $BAMCHECK == 1 ]; then
    for f in $FASTQFILES; do
        FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
        SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
        awk -v b="BAMCHECK=$BAMCHECK;" '{ if(NR==23) {print b} else {print $0} } ' $SHSCRIPT > "$SHSCRIPT.temp"; mv "$SHSCRIPT.temp" $SHSCRIPT;
    done
else continue;
fi

#CUFFLINKS
while true; do
    read -p "Do you want to pipe directly to cufflinks after completion (y or n)? " yn
        case $yn in
        [Yy]* ) CUFFLCHECK=1; break;;
        [Nn]* ) CUFFLCHECK=0; break;;
        * ) echo "\n"; echo "Please answer y or n: ";
    esac
done
#WRITE VARIABLES TO CUFFLINK.sh
WAIT=/users/j/a/jargordo/bin/scripts/wait_for_accepted_hits.sh
if [ $CUFFLCHECK == 1 ]; then
printf "\nWriting CUFFLINKS script... "
    for f in $FASTQFILES; do
        FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
        AHITSDIR=$PATHTOWD/$FNAME
        CUFFLSCRIPT="submit_CUFFLINKS_""$FNAME"".sh"; cp /users/j/a/jargordo/bin/scripts/run_CUFFLINKS_on_VACC.sh $CUFFLSCRIPT
        WAIT="wait_for_$FNAME.sh"
        JOBNAME="#PBS -N ""$FNAME""_CUFFLINKS"; awk -v a="$JOBNAME" '{ if(NR==7) { print a} else {print $0} } ' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v b="WD=$PATHTOWD;" '{ if(NR==8) {print b} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v c="INPUT=$AHITSDIR;" '{ if(NR==9) {print c} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v d="SPECIES=$SPECIES;" '{ if(NR==10) {print d} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v e="REFERENCE=$REFERENCE;" '{ if(NR==11) {print e} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v f="LIBTYPE=$LIBTYPE;" '{ if(NR==12) {print f} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
        awk -v g="WAIT=$WAIT;" '{ if(NR==13) {print g} else {print $0} }' $CUFFLSCRIPT > "$CUFFLSCRIPT.temp"; mv "$CUFFLSCRIPT.temp" $CUFFLSCRIPT
#wait scripts
        awk -v a="FNAME=$FNAME;" -v b="WD=$WD;" -v c="CUFFLSCRIPT=$CUFFLSCRIPT;" '{ if(NR==4) { print a,b,c } else {print $0} } ' /users/j/a/jargordo/bin/scripts/wait_for_accepted_hits.sh > $WAIT
        chmod a+x $WAIT
        nohup sh $WAIT >> temp.out &
    done
sleep 1
printf "done. \n"
fi

#HTSEQ
#CUFFDIFF
#printf "WARNING: Things get trickier here because there are multiple conditions/files in play. \n"
#printf "If you select CUFFDIFF you have to be very carefull about file names and conditions for it to work. \n"
#printf "\n"
#while true; do
#    read -p "Do you want to pipe results directly to CUFDIFF after completion: (y or n)? " yn
#        case $yn in
#        [Yy]* ) DIFFCHECK=1; break;;
#        [Nn]* ) DIFFCHECK=0; break;;
#        * ) echo "\n"; echo "Please answer yes or no.";
#    esac
#done
#if [ $DIFFCHECK == 1 ]; then
#fi

printf "\n";
while true; do
    read -p "Last chance to bail (q) or submit job (s)? " qs
        case $qs in
        [Qq]* ) printf "\nScripts were written but not submitted (qsub)... exiting"; exit; break;;
        [Ss]* ) break;;
        * ) echo "\n"; echo "Please select quit (q), print options (o) or submit (s): ";
    esac
done

printf "\n";
printf "No more questions. \n";
sleep 1
printf "Dramatic pause.... \n";
printf "\n";
sleep 1
printf "submitting jobs  \n";
printf "\n";

#FINAL QSUB
for f in $FASTQFILES; do
    FNAME=$(echo $f | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
    SHSCRIPT="submit_TOPHAT_""$FNAME"".sh"
    chmod a+x $SHSCRIPT
    echo `qsub $SHSCRIPT`;
done