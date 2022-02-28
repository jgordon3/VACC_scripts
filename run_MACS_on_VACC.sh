#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=24:00:00
#PBS -N MACS
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y_%H%M%S`
WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
mkdir $WD
cd $WD

#---------SET_CHROM_FOR_WIGTOBIGWIG---------
CHROMSIZE="$HOME/scratch/indexes/hg19.chrom.sizes"

#----------FILES_TO_PROCESS-----------------------
PARENTFOLDER=/users/j/a/jargordo/scratch/BOWTIE_06_06_2013_154325/
INPUT1=MCF-7_input_output.bowtie
INPUT2=MCF-10a_input_output.bowtie
INPUT3=MDA_input_output.bowtie

AA=MCF7H3K9ME2_R1.output.bowtie
BB=MCF7H3K27AC_R1.output.bowtie
CC=MCF7H4K20ME3_R1.output.bowtie
DD=MCF7H3K27ME3_R1.output.bowtie
EE=MCF7H3K4AC_.output.bowtie

FF=MCF10AH3K27AC_combined.output.bowtie
GG=MCF10AH3K27ME3_combined.output.bowtie
HH=MCF10AH4K20ME3_combined.output.bowtie

II=MDAH3K27ME3_R1.output.bowtie
JJ=MDAH4K20ME3_R1.output.bowtie
KK=MDAH3K27AC_R1.output.bowtie

#---------------LOOPS----------------------------
for sample in $AA; do
BOWTIE_TESTFILE="$PARENTFOLDER""$sample"
BOWTIE_CONTROLFILE="$PARENTFOLDER""$INPUT1"
ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
echo "$ABRV_NAME"
echo "------------------------MACS-----------------------------------"
macs -t $BOWTIE_TESTFILE -c $BOWTIE_CONTROLFILE  -f BOWTIE --name=$ABRV_NAME --bw=400 --pvalue=1.00e-5 -w --single-profile
done

for sample in $FF $GG $HH; do
BOWTIE_TESTFILE="$PARENTFOLDER""$sample"
BOWTIE_CONTROLFILE="$PARENTFOLDER""$INPUT2"
ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
echo "$ABRV_NAME"
echo "------------------------MACS-----------------------------------"
macs -t $BOWTIE_TESTFILE -c $BOWTIE_CONTROLFILE  -f BOWTIE --name=$ABRV_NAME --bw=400 --pvalue=1.00e-5 -w --single-profile
done

for sample in $II $JJ $KK; do
BOWTIE_TESTFILE="$PARENTFOLDER""$sample"
BOWTIE_CONTROLFILE="$PARENTFOLDER""$INPUT3"
ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
echo "$ABRV_NAME"
echo "------------------------MACS-----------------------------------"
macs -t $BOWTIE_TESTFILE -c $BOWTIE_CONTROLFILE  -f BOWTIE --name=$ABRV_NAME --bw=400 --pvalue=1.00e-5 -w --single-profile
done

for  sample in $AA $BB $CC $DD $EE $FF $GG $HH $II $JJ $KK; do
ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
BTBB_INFILE="$ABRV_NAME""_peaks.bed"
BTBB_OUTFILE="$ABRV_NAME""_peaks.bb"
WTBW_DIR="./""$ABRV_NAME""_MACS_wiggle/"
WTBW_INFILE="$WTBW_DIR""treat/""$ABRV_NAME""_treat_afterfiting_all.wig"
WTBW_OUTFILE="$ABRV_NAME""_signal.bw"
echo "$ABRV_NAME"
echo "------------------------WIGTOBIGWIG/BEDTOBIGBED-----------------------------------"
gunzip "$WTBW_INFILE"".gz"
wigToBigWig -clip $WTBW_INFILE $CHROMSIZE $WTBW_OUTFILE
bedToBigBed -type=bed3+2 $BTBB_INFILE $CHROMSIZE $BTBB_OUTFILE
rm -r $WTBW_DIR
echo "track type=bigWig" "name=$ABRV_NAME""_sig" "description="'"'"MACS big wig from $ABRV_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/$WTBW_OUTFILE" "autoScale=off viewLimits=1:20 smoothingWindow=2 windowingFunction=maximum" "color=$((RANDOM %256 + 1)),$((RANDOM %256 +1)),$((RANDOM %256 +1))" "visibility=dense" >> UCSC_CONTROL_FILE_$DATE.txt
echo "track type=bigBed" "name=$ABRV_NAME""_peaks" "description="'"'"peaks from $ABRV_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/$BTBB_OUTFILE" >> UCSC_CONTROL_FILE_$DATE.txt
done

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
which macs

# MACS2

#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID