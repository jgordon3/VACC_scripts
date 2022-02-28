#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=4gb,pvmem=5gb
#PBS -l walltime=24:00:00
#PBS -N FIMO_CHROM_SPLIT
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

#Script splits chromosomes into small parts (100Kb) to run FIMO motifs across genome and assembles a bed file as output

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
#WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE"_"$PBS_O_JOBID
WD=/users/j/a/jargordo/scratch/FIMO_TEST/
if [ -d "$WD" ]
then
cd $WD
else
mkdir $WD
cd $WD
fi

GENOME="mm10"
#set PVALUE, MOTIF, etc
PVALUE=1e-3
MOTIF=Runx2_motif.txt
MOTIFNAME=$(echo "$f2" | awk -F '.' '{print $1}')
OUTPUT_FILE=$GENOME"_"$MOTIFNAME"_"$PVALUE".bed"
CHROMDIR="/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes"
cp $CHROMDIR/* ./
echo "copying chromosomes from $CHROMDIR"
# loop for every chromosome
CHROM="$WD""*.fa"

for f in $CHROM; do
#count lines in file
LINES=$(wc -l < $f)
NFILES=$(($LINES/100000))
FNAME=$(echo "$f" | awk -F '/' '{print $13}')

#split fasta into 100K segments using pyfasta, seprate into NFILES number of files
pyfasta split -n $NFILES -k 100000 $f
echo "$f split finished"
#make list of chrom split files
CHROM_SPLIT="$WD""*100Kmer*"

# loop 2
for f2 in $CHROM_SPLIT; do
#extract the file interation
MPLYER=$(echo "$f2" | awk -F '.' '{print $2}')
if [ "$MPLYER" = "00" ]
then
MPLYER="0"
else
MPLYER=`echo $MPLYER |sed 's/^0*//'`
fi

#add 100000K * MPLYER variable to each bed record
ADJUST=$(($MPLYER*100000))
#run FIMO
fimo --thresh $PVALUE $MOTIF $f2
#create bed from FIMO txt file
tail -n +2 ./fimo_out/fimo.txt > holder1.txtx
awk '{split($2,chr,"_");print chr[1]"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' holder1.txtx > holder2.txtx

#adjust chromosome values on bed
awk -v ad=$ADJUST '{print $1"\t"($2 + ad)"\t"($3 + ad)"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' holder2.txtx > adjusted_fimo.txtx
sort -k1,1 -k2,2n adjusted_fimo.txtx >> $OUTPUT_FILE
#clean up
rm *.txtx
rm -r fimo_out
echo "done run $MPLYER"
#end loop 2
done
#clean up
rm $CHROM_SPLIT
rm *.flat
rm *.gdx

done

#clean up splits and beds

