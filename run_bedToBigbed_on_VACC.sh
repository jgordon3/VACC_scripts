#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=24:00:00
#PBS -N BEDTOBIGBED
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

#REQUIRES: SPECIES, INPUT,

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
#WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE"_"$PBS_O_JOBID
if [ -d $WD ]; then cd $WD; else mkdir $WD; cd $WD; fi

#---------SET_CHROM_FOR_WIGTOBIGWIG---------
if [ -z $SPECIES ];
then SPECIES="mm10"; echo "GENOME HAS BEEN SET TO $SPECIES BY DEFAULT"
else echo "GENOME WAS SET TO $SPECIES"
fi

if  [ $SPECIES = "mm10" ]; then
CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/mm10.chrom.sizes"
elif [ $SPECIES = "hg19" ]; then
CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/GenomeStudio/Homo_sapiens/UCSC-hg19/hg19.chrom.sizes"
else echo "$GENOME is not set";
fi

FILE=$(echo "$INPUT" | awk -F '/' '{print $NF}')
NAME=$(echo "$FILE" | awk -F '.' '{print $1}')
SORTED="$PBS_JOBID"".sorted"
CLIPPED="$PBS_JOBID"".clipped"
BIGBED="$NAME"".bb"
FBED="$NAME""_sorted.bed"
sort -k1,1 -k2,2n $FILE > $SORTED
bedClip -verbose=2 $SORTED $CHROMSIZE $CLIPPED
bedToBigBed -type=bed3+2 $CLIPPED $CHROMSIZE $BIGBED
mv $CLIPPED $FBED
rm $SORTED

echo $FILE
echo $BIGBED
echo $FBED

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`



#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID