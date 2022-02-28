#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=24:00:00
#PBS -N WIGTOBIGWIG
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

#---------SET_CHROM_FOR_WIGTOBIGWIG---------
CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/GenomeStudio/Homo_sapiens/UCSC-hg19/hg19.chrom.sizes"

#----------FILES_TO_PROCESS-----------------------
PARENTFOLDER=/users/j/a/jargordo/scratch/Dan_zlab/
cd $PARENTFOLDER
FILES=/users/j/a/jargordo/scratch/Dan_zlab/*.wig

#AA=Encode_MCF7_H3K27me3_R1.bowtie
#BB=Encode_MCF7_H3K27me3_R2.bowtie
#CC=Encode_MCF7_H3K4me3_R1.bowtie
#DD=Encode_MCF7_H3K4me3_R2.bowtie
#HH=MCF-7_h3k4me3_raw.bowtie

for sample in $FILES; do
FILENAME=$(echo "$FILES" | awk -F '/' '{print $8}')
ABRV_NAME=$(echo "$FILENAME" | awk -F '.' '{print $1}')
echo "$ABRV_NAME"

# edit this if changed
MACS2_BEDGRAPH="$ABRV_NAME""_treat_pileup.bdg"
MACS2_CONTROL="$ABRV_NAME""_control_lambda.bdg"
MACS2_OUTFILE1="$ABRV_NAME""_fold_enrichment.bdg"
MACS2_OUTFILE2="$ABRV_NAME""_logLR.bdg"

GRAPHCOLOUR=$((RANDOM %256 + 1)),$((RANDOM %256 +1)),$((RANDOM %256 +1))

for outfile in $MACS2_BEDGRAPH $MACS2_CONTROL $MACS2_OUTFILE1 $MACS2_OUTFILE2; do
FINAL_NAME=$(echo "$outfile" | awk -F '.' '{print $1}')
CLIPPED="$FINAL_NAME"".clipped"
FINALWIG="$FINAL_NAME"".bw"
TARNAME="$outfile"".tar.gz"
echo "---------------bedGraphToBigWig_for_$outfie-----------------------------------"
bedClip -verbose=2 $outfile $CHROMSIZE $CLIPPED
bedGraphToBigWig $CLIPPED $CHROMSIZE $FINALWIG
rm *.clipped

#OPTION: remove bedgraphs?
#rm $outfile
#OPTION: archive edgraphs?
#tar -zcvf $TARNAME $outfile

#7A)-----GENERATE_CONTROL_FILE_FOR_UCSC_GENOME_BROWSER-------------------------

echo "track type=bigWig" "name=$FINAL_NAME" "description="'"'"MACS2 big wig from $FINAL_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/$FINALWIG" "autoScale=off viewLimits=1:20 smoothingWindow=2 windowingFunction=maximum" "maxHeightPixels=128:48:8" "color=$GRAPHCOLOUR" "visibility=dense" >> UCSC_CONTROL_FILE_BIGWIGS_$DATE.txt
done

#8)-----CONVERT_BEDS_TO_BIGBED-------------------------

PEAKSBED="$ABRV_NAME""_peaks.bed"
SUMMITBED="$ABRV_NAME""_summits.bed"
for outfile1 in $PEAKSBED $SUMMITBED; do
FINAL_NAME=$(echo "$outfile1" | awk -F '.' '{print $1}')
CLIPPEDBED="$FINAL_NAME"".clipped"
FINALBED="$FINAL_NAME"".bb"
bedClip -verbose=2 $outfile1 $CHROMSIZE $CLIPPEDBED
bedGraphToBigWig $CLIPPEDBED $CHROMSIZE $FINALBED

#8A)-----GENERATE_CONTROL_FILE_FOR_UCSC_GENOME_BROWSER-------------------------
echo "track type=bigBed" "name=$FINAL_NAME" "description="'"'"peaks/summits from $FINAL_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/$FINALBED" >> UCSC_CONTROL_FILE_BIGBED_$DATE.txt
done

done

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
#which fastqc
#which bowtie


#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID