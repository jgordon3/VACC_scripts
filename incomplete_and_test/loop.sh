#!/bin/bash
while read AA BB CC DD
do

JOB=`qsub -m bea -N TESTRUN - << EOJ

#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=1gb,pvmem=1gb
#PBS -l walltime=03:00:00
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu

WD=$HOME/scratch/TESTRUN
if [ -d "$WD" ]
then
cd $WD
else
mkdir $WD
cd $WD
fi

#----------BOWTIE_INDEX_--------------------------
#INDEXPATH=/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
INDEXPATH=/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome

echo ${AA} ${BB} ${CC} ${DD}
echo $INDEXPATH
EOJ
`

echo "JobID = ${JOB} for parameters ${AA} ${BB} ${CC} ${DD} submitted on `date`"
done < configure.txt
exit
