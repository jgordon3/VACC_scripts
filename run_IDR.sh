#!/bin/sh

# runs IDR pipeline for ChIP-seq replicate analysis
# Packages can be found here:
# http://cran.r-project.org/web/packages/idr/index.html
# https://sites.google.com/site/anshulkundaje/projects/idr
# REQUIRES: genome_table.txt (same as ChromInfo.txt) ($CHROMSIZE)
# REQUIRES: peakfile 1, peakfile 2

#FOLDER CONTAINING IDR CODE
IDR="$HOME/Desktop/idrCode/*"

WD="$HOME/Desktop/RUNX2_IDR_ANALYSIS_04_16_2014/"
if [ -d $WD ]; then cd $WD; else mkdir $WD; cd $WD; fi
cp $IDR $WD



FOLDER=/Volumes/SteinLab-Research\$/Jonathan/transcription_factor_paper_data/BMSC_RUNX2_MACS_04_15_2014/
FOLDERPEAKS="$FOLDER/BMSC*_p_1E-03_peaks.narrowPeak"

# USE MACS TO GENERATE 1e-03 PVALUE PEAKS
for f in $FOLDERPEAKS; do
FILE1=$f
FILEN=$(echo "$f" | awk -F '/' '{print $NF}')
SFILEN=$(echo "$FILEN" | awk -F '_' '{print $1 "_" $2}')
DIR1="$WD""$SFILEN"
mkdir "$DIR1"
ODIR="$DIR1/""$SFILEN"

R1NAME="$SFILEN""_R1"
R2NAME="$SFILEN""_R2"
VERSUS="$ODIR""$R1NAME""vs""$R2NAME"
FILET=$"$FOLDER2/""$R2NAME*_p_1E-03_peaks.narrowPeak"
FILE2=$(echo $FILET)
scp $VACC/scratch/indexes/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/ChromInfo.txt $WD/genome_table.txt
printf "$R1NAME $R2NAME \n"
Rscript batch-consistency-analysis.r $FILE1 $FILE2 -1 $ODIR 0 F p.value
#Rscript batch-consistency-plot.r [npairs] [output.prefix] [input.file.prefix1] [input.file.prefix2] [input.file.prefix3]
Rscript batch-consistency-plot.r 1 $VERSUS $ODIR
TOTAL_PEAKS1=$(cat $FILE1 | wc -l); echo " TOTAL PEAKS FILE 1: $TOTAL_PEAKS1"
TOTAL_PEAKS2=$(cat $FILE2 | wc -l); echo " TOTAL PEAKS FILE 2: $TOTAL_PEAKS2"
OLAPS="$ODIR-overlapped-peaks.txt"
OVER_PEAKS=$(cat $OLAPS | awk '$11 <= 0.01 {print $0}' | wc -l);echo " OVERLAPPED PEAKS: $OVER_PEAKS"
done
