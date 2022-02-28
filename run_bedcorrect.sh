#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=24:00:00
#PBS -N BedCorrect
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
#WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
WD=/users/j/a/jargordo/scratch/Dan_zlab
if [ -d "$WD" ]
then
cd $WD
else
mkdir $WD
cd $WD
fi

mkdir Corr

for sample in SteinAML1ETO_CS_Comb_Peaks.bed SteinAML1ETO_Comb_1and3_Peaks.bed SteinAML1ETO_Comb_2and3_Peaks.bed SteinAML1ETO_Comb_Peaks.bed SteinAML1ETO_Rep1_Peaks.bed SteinAML1ETO_CS_Rep1_Peaks.bed SteinAML1ETO_CS_Rep2_Peaks.bed SteinAML1ETO_Rep2_Peaks.bed SteinAML1ETO_insilico_Peaks.bed SteinAML1_AB_Rep1_Peaks.bed SteinAML1_consensus.bed SteinAML1_AB_Comb_Peaks.bed SteinAML1_AB_Rep2_Peaks.bed SteinAML1_AM_Comb_Peaks.bed SteinAML1_AM_Rep2_Peaks.bed SteinAML1_insilico_Comb_Peaks.bed SteinH3K27me3_Comb_Peaks.bed SteinH3K27me3_Rep1_Peaks.bed SteinH3K27me3_Rep2_Peaks.bed SteinH3K4me3_Comb_Peaks.bed SteinH3K4me3_Rep1_Peaks.bed SteinH3K4me3_Rep2_Peaks.bed SteinH3K9ac_Rep1_Peaks.bed SteinNcor_Comb_Peaks.bed SteinNcor_Rep1_Peaks.bed SteinNcor_Rep2_Peaks.bed Steinp300_Comb_Peaks.bed Steinp300_Comb_Rep1and3_Peaks.bed Steinp300_Rep1_Peaks.bed Steinp300_Rep2_Peaks.bed Steinp300_Rep3_Peaks.bed siH3K9ac_Rep1_Peaks.bed SteinAML1ETO_Rep3_Peaks.bed SteinAML1_PWM.bed SteinAML1_AM_Rep1_Peaks.bed; do

ABRV_NAME=$(echo "$sample" | awk -F '.' '{print $1}')
OUTPUT_NAME=$ABRV_NAME"_corr.bed"

cut -f 2-8 $sample| sed 1d | sort -k1,1 -k2,2n > ./Corr/$OUTPUT_NAME

done

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
which bowtie

#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID