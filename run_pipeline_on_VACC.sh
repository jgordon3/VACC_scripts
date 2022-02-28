#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=24:00:00
#PBS -N BOWTIE
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
#WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
WD=$HOME/scratch/BOWTIE_07_25_2013
if [ -d "$WD" ]
then
cd $WD
else
mkdir $WD
cd $WD
fi

#----------FILES_TO_PROCESS-----------------------
PARENTDIR=/users/j/a/jargordo/scratch/sequences/H3K4ME3_comp_07_23_2013/

AA=Encode_MCF7_H3K27me3_R1.fastq
BB=Encode_MCF7_H3K27me3_R2.fastq
#CC=Encode_MCF7_H3K4me3_R1.fastq
#DD=Encode_MCF7_H3K4me3_R2.fastq
#EE=Encode_MCF7_Input_R1.fastq
#FF=Encode_MCF7_Input_R1UW.fastq
#GG=Encode_MCF7_Input_R2.fastq
HH=MCF-7_h3k4me3_raw.fastq
II=MCF-7_input_raw.fastq

#----------BOWTIE_INDEX_--------------------------
echo "--------------Index_used:UCSC_hg19_2013-03-06-11-23-03_for_Bowtie2------------------"
INDEXPATH=/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome

for sample in $HH $II; do


#_____gunzip?______
gunzip "$PARENTDIR$sample"

INPUT_NAME="$PARENTDIR""$sample"
ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
OUTPUT_NAME=$ABRV_NAME".bowtie"

echo "------------------------BOWTIE_for_$ABRV_NAME-----------------------------------"
bowtie2 -N 0 -L 28 -x $INDEXPATH -U $INPUT_NAME -S $OUTPUT_NAME

#gzip "$INPUT_NAME"
done

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
which bowtie

#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID