#!/bin/bash
#PBS -l nodes=4:ppn=8,pmem=8gb,pvmem=9gb
#PBS -l walltime=03:00:00
#PBS -N ChIP_SEQ_PIPELINE
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y_%H%M%S`
WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
mkdir -p $WD/{fastqc,bowtie,macs,logs}
cd $WD

#options: hg19, mm9, mm10
INDEXPATH="$HOME/bin/bowtie-0.12.9/indexes/mm9"

# edit this if changed
PATHTOSEQFILE="$HOME/scratch/sequences/s_4_sequence.txt"

#echo "------------------------FASTQC-----------------------------------"
#fastqc $PATHTOSEQFILE --noextract -o $WD/fastqc


# test paralell 
echo "------------------------BOWTIE-----------------------------------"

cat $PATHTOSEQFILE | parallel --block 300k --recstart '>' --pipe bowtie -n 2 -l 36 $INDEXPATH >> $WD/bowtie/test_output.bowtie

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
#which fastqc
which bowtie
echo "sequence" $PATHTOSEQFILE
echo $PBS_O_JOBID
mv $HOME/bin/scripts/"$PBS_JOBNAME".o"$PBS_O_JOBID" $WD/logs/

#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID