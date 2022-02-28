#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=04:00:00
#PBS -N BOWTIE
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

# DO NOT MODIFY UNLESS YOU KNOW WHAt YOU ARE DOING
# Use submit script (submit_bowtie.sh) to submit this script or try "qsub -v INPUT=/path/to/fastq,GENOME=mm10 run_bowtie_on_VACC.sh" to do it directly
# REQUIRES: INPUT (A input fastq)
# REQUIRES: GENOME a reference genome (e.g. mm10 or hg19)
# OPTION: IFN: Include folder name (Y or N). Include the sub folder in the naming convention (e.g. H3K4me3/replicate_file.fastq would be renamed to H3K4me3_replicate_file.fastq)
# OPTION: FNC: File naming convention (Y or N): File is named in CELL_CONDITION_REPLICATE format. This will use
# OPTION: TRIM: Trim 5' and 3' ends. Defalt is -5 8 -3 35 (effective length is 57 bp)
# OPTION: MISMATCH
# OPTION: LENGTH
# OPTION: FILTER

# don't change this unless absolutely needed
DATE=`date +%m_%d_%Y`
#WD=/users/j/a/jargordo/scratch/ENCODE_alignments_11_11_2013/
if [ -d "$WD" ]; then cd $WD; else mkdir WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE; cd $WD; fi

#----------SET_BOWTIE_INDEX_--------------------------
# from the submit script
if [ $GENOME = "mm10" ]; then INDEXPATH=/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/Sequence/Bowtie2Index/genome;
elif [ $GENOME = "hg19" ]; then INDEXPATH=/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome;
else echo "GENOME could not be set"; exit 1;
fi

#----------FILES_TO_PROCESS-----------------------
#CHECK FOR GZIP
GZIP=$(echo $INPUT | awk -F '.' '{print $NF}')
if [ $GZIP = "gz" ]; then gunzip $INPUT; INPUT=$(echo $INPUT | rev | cut -d. -f2- | rev); else echo "not gzipped"; fi

#DEFINE INPUTS AND OUTPUTS
FILE=$(echo $INPUT | awk -F '/' '{print $NF}')
ABRV_NAME=$(echo "$FILE" | awk -F '.' '{print $1}')

#CHECK IF FNC AND IFN ARE SET AND NAME APPROPRIATELY
if [ $FNC = "Y" ]; then
    echo "FILE NAME CONVENTION SET."
    CELL=$(echo "$ABRV_NAME" | awk -F '_' '{print $1}')
    echo "CELL: $CELL"
    CONDITION=$(echo "$ABRV_NAME" | awk -F '_' '{print $2}')
    echo "CONDITION: $CONDITION"
    REPLICATE=$(echo "$ABRV_NAME" | awk -F '_' '{print $3}')
    echo "REPLICATE: $REPLICATE"
    ABRV_NAME=$CELL"_"$CONDITION"_"$REPLICATE
    echo $ABRV_NAME
else echo "FILE NAME CONVENTION NOT SET. YOU ARE STUCK WITH THE LONG FILE NAME"; fi

#INCLUDE FOLDER NAME
#if [ $IFN = "Y" ]; then FOLNAME=$(echo $INPUT| awk -F '/' '{print $(NF-1)}'); OUTPUT_NAME="$FOLNAME""_""$ABRV_NAME"; else OUTPUT_NAME="$ABRV_NAME"; fi

# OPTIONAL PREFILTERING USING TRIMMOMATIC
if [ FILTER = "Y" ]; then
FFASTQ="$ABRV_NAME""_filtered.fastq";
java -jar /users/j/a/jargordo/bin/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 $INPUT $FFASTQ ILLUMINACLIP:/users/j/a/jargordo/bin/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36;
INPUT=$FFASTQ
else echo "NO PREFILTER"
fi

# CHECK IF OPTIONAL BOWTIE PARAMETERS ARE SET, IF NOT USE DEFAULTS:
if [ -z $TRIM5 ]; then TRIM5="8"; else echo ""; fi
if [ -z $TRIM3 ]; then TRIM3="20"; else echo ""; fi
if [ -z $MISMATCH ]; then MISMATCH="0"; else echo ""; fi
if [ -z $LENGTH ]; then LENGTH="22"; else echo ""; fi

echo "Input: $INPUT"
echo "Output: $OUTPUT_NAME"
echo "Index: $INDEXPATH "

# ACTUAL BOWTIE CALLS HAPPENS HERE
bowtie2 -5 $TRIM5 -3 $TRIM3 -N $MISMATCH -L $LENGTH -x $INDEXPATH -U $INPUT -S $OUTPUT_NAME
echo "Bowtie call: bowtie2 -5 $TRIM5 -3 $TRIM3 -N $MISMATCH -L $LENGTH -x $INDEXPATH -U $INPUT_NAME -S $OUTPUT_NAME"


#MOVING STUFF AROUND AND RENAMING
# count alignments and add to file name
ALIGNS=$(wc -l "$OUTPUT_NAME"| awk '{print $1}')
BOWTIE_OUTPUT=$OUTPUT_NAME"_"$ALIGNS".bowtie"

#name BAM and TAGALIGN files
BAMFILE=$OUTPUT_NAME"_"$ALIGNS".bam"
SORTEDBAM=$OUTPUT_NAME"_"$ALIGNS".sorted"
TAGALIGN=$OUTPUT_NAME"_"$ALIGNS".tagAlign.gz"

mv $OUTPUT_NAME $BOWTIE_OUTPUT

#gzip $INPUT_NAME
#gzip $BOWTIE_OUTPUT

#CONVERT SAM TO BAM AND TAGALIGN FOR IDR
# quality filter (-q) and generate a BAM file
samtools view -bS -o $BAMFILE $BOWTIE_OUTPUT
#write tagAlign
bamToBed -i $BAMFILE | awk 'BEGIN{FS="\t";OFS="\t"}{$4="N"; print $0}' | gzip -c > $TAGALIGN

#Calculate unique alignments (MAPQ=1 is unique) and quality (MAPQ > 30)
UNIQUES=$(samtools view -Sq 1 $BOWTIE_OUTPUT | wc -l)
QUALITY=$(samtools view -Sq 30 $BOWTIE_OUTPUT | wc -l)
echo "reads with unique alignments (from the BAM file): $UNIQUES"
echo "reads with MAPQ quality scores above 30 (from the BAM file): $QUALITY"
# replace BAM file with a sorted BAM file
samtools sort $BAMFILE $SORTEDBAM

#PICARD TOOLS
#$PICARDPATH="/users/j/a/jargordo/bin/picard-tools-1.107/"
java -jar /users/j/a/jargordo/bin/picard-tools-1.107/QualityScoreDistribution.jar INPUT=$BAMFILE OUTPUT=$OUTPUT_NAME"_Quality_Dist" CHART_OUTPUT=$OUTPUT_NAME"_Quality_Dist.pdf"
rm $BAMFILE
echo "Bowtie output (for MACS2):$BOWTIE_OUTPUT"
echo "TagAlign file (for SPP, IDR): $TAGALIGN"
echo "BAM file (for whatever): $SORTEDBAM"

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`
which bowtie2


#Here's the command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID