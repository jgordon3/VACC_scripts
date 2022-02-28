#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=30:00:00
#PBS -N FASTX
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
if [ -d "$WD" ]
then
cd $WD
else
mkdir $WD
cd $WD
fi

#----------FILES_TO_PROCESS-----------------------
PARENTFOLDER=/users/j/a/jargordo/scratch/BOWTIE_REALIGNMENTS_09_24_2013/

AA=MCF10A_H3K4me3_ATGTCA_L007_CONTROL_R1_000.fastq

for sample in $AA; do

INPUT="$PARENTFOLDER$sample"
OUTPUT=$(echo "$sample" | awk -F '.' '{print $1 "_collapsed.fastq" }')

fastx_collapser -v -Q 33 -i $INPUT -o $OUTPUT

done
