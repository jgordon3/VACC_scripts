#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea
#REQUIRES: SPECIES
#OPTIONS:WD (working directory), REFERENCE (will be set to NCBI by default), LIBTYPE (Library type options for TopHat and Cufflinks: default is unstranded)
# SEE: Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., … Pachter, L. (2012).
# Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.
# Nature Protocols, 7(3), 562–78. doi:10.1038/nprot.2012.016
# don't change this script unless absolutely needed
#THIS SPACE IS BLANK FOR ADDITIONAL VARIABLES

#----------------------------
# SETTING DIRECTORIES
DATE=`date +%m_%d_%Y`
if [ -z $WD ]; then WD=$HOME/scratch/"TOPHAT_"$DATE; else printf "."; fi
if [ -d "$WD" ]; then cd $WD; else mkdir $WD; cd $WD; fi

#---------LOCATION_OF_GENE_REFERENCE_FILES---------
#need to finish
if [ -z $SPECIES ]; then SPECIES="mm10"; echo "GENOME HAS BEEN SET TO $SPECIES BY DEFAULT"; else printf "."; fi
if [ -z $REFERENCE ]; then REFERENCE="UCSC"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else printf "."; fi
if [ $SPECIES = "mm10" ]; then FULLSPECIES="Mus_musculus"; elif [ $SPECIES = "hg19" ]; then FULLSPECIES="Homo_sapiens"; else echo "REFERENCE: $REFERENCE COULD NOT BE SET"; fi

INDEXFOLDER=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES

if [ -z $REFERENCE ]; then REFERENCE="UCSC"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else printf "."; fi
    if [ $REFERENCE = "UCSC" ]; then
        GENES=$INDEXFOLDER/UCSC/$SPECIES/Annotation/Genes/genes.gtf
        GENOME=$INDEXFOLDER/UCSC/$SPECIES/Sequence/Bowtie2Index/genome
        GENOMEFA=$INDEXFOLDER/UCSC/$SPECIES/Sequence/Bowtie2Index/genome.fa
        CHROMSIZE=$INDEXFOLDER/UCSC/$SPECIES/Annotation/Genes/ChromInfo.txt
    elif [ $REFERENCE = "NCBI" ]; then
        GENES=$INDEXFOLDER/NCBI/Annotation/Genes/genes.gtf
        GENOME=$INDEXFOLDER/NCBI/Sequence/Bowtie2Index/genome
        GENOMEFA=$INDEXFOLDER/NCBI/Sequence/WholeGenomeFasta/genome.fa
        CHROMSIZE=$INDEXFOLDER/NCBI/GenomeStudio/$FULLSPECIES/ChromInfo.txt
    elif [ $REFERENCE = "ENSEMBL" ]; then
        GENES=$INDEXFOLDER/ENSEMBL/Annotation/Genes/genes.gtf
        GENOME=$INDEXFOLDER/ENSEMBL/Sequence/Bowtie2Index/genome
        GENOMEFA=$INDEXFOLDER/ENSEMBL/Sequence/WholeGenomeFasta/genome.fa
        CHROMSIZE=$INDEXFOLDER/GenomeStudio/$FULLSPECIES/ChromInfo.txt
    elif [ $REFERENCE = "GENCODE" ]; then
        GENES=$INDEXFOLDER/GENCODE/Annotation/Genes/genes.gtf
        GENOME=$INDEXFOLDER/UCSC/$SPECIES/Sequence/Bowtie2Index/genome
        GENOMEFA=$INDEXFOLDER/UCSC/$SPECIES/Sequence/Bowtie2Index/genome.fa
        CHROMSIZE=$INDEXFOLDER/UCSC/$SPECIES/Annotation/Genes/ChromInfo.txt
    else echo "REFERENCE: $REFERENCE COULD NOT BE SET"
fi

# DEFINE INPUTS/OUTPUTS
printf "\n"
echo "PATHTOFASTQ: $PATHTOFASTQ"
# CHECK FOR GZIP
GZIP=$(echo $PATHTOFASTQ | awk -F '.' '{print $NF}')
if [ $GZIP = "gz" ]; then gunzip $PATHTOFASTQ; PATHTOFASTQ=$(echo $PATHTOFASTQ| rev | cut -d. -f2- | rev); else printf "."; fi
while  [  ! -f $PATHTOFASTQ ]; do sleep 120; done

#FILE NAMING STUFF
#inputs
FILE=$(echo $PATHTOFASTQ | awk -F '/' '{print $NF}')
ABRV_NAME=$(echo "$FILE" | awk -F '.' '{print toupper($1)}')
#outputs

#TRIM READS
if [ ! -z "$TRIMCHECK" ]; then
INPUT=$PATHTOFASTQ;
OUTPUT="$ABRV_NAME""_trimmed.fastq"
java -jar /users/j/a/jargordo/bin/Trimmomatic-0.32/trimmomatic-0.32.jar SE -phred33 $INPUT $OUTPUT ILLUMINACLIP:/users/j/a/jargordo/bin/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:$LEADING TRAILING:$TRAILING SLIDINGWINDOW:$SLIDINGWINDOW MINLEN:36; sleep 900;
PATHTOFASTQ=$OUTPUT
fi

#BUILD TOPHAT EXPRESSION
TOPHATPARAM="tophat -G $GENES -o $ABRV_NAME"

# CHECK FOR LIBTYPE
if [ -z $LIBTYPE ]; then LIBTYPE="fr-unstranded"; fi
TOPHATPARAM="$TOPHATPARAM --library-type=$LIBTYPE"
# CHECK FOR ALENGTH
if [ ! -z $ALENGTH ]; then echo "Anchor length was changed to $ALENGTH (default: 8)"; TOPHATPARAM="$TOPHATPARAM -a $ALENGTH"; fi
# CHECK FOR JUNCTIONS
if [ ! -z $JUNCTIONS ]; then echo "Junctions were changed to: $JUNCTIONS"; TOPHATPARAM="$TOPHATPARAM $JUNCTIONS"; fi
# CHECK FOR JUNCTION FILE
if [ ! -z $JUNCTIONF ]; then echo "Junctions file: $JUNCTIONF was used to assign junctions"; TOPHATPARAM="$TOPHATPARAM -j $JUNCTIONF"; fi
#eval TOPHATPARAM=$TOPHATPARAM
TOPHATPARAM="$TOPHATPARAM $GENOME";
TOPHATPARAM="$TOPHATPARAM $PATHTOFASTQ"



echo "WD: $WD"
echo "GENOME: $SPECIES"
if [ ! -z "OUTPUT" ]; then echo "FASTQ WAS TRIMMED USING TRIMMOMATIC"; fi
printf "\n"
echo "TOPHAT SUBMITED WITH FOLLOWING PARAMETERS: \n"
echo "$TOPHATPARAM"

#this is the submit line
echo `$TOPHATPARAM`


if [ $BAMCHECK == 1 ]; then
#SORT --> BAM --> BEGRAPH --> BIGWIG
CUFFIN=$ABRV_NAME"/accepted_hits.bam"
SBAM=$ABRV_NAME"/accepted_hits_sorted"
PLUS=$ABRV_NAME/$ABRV_NAME"_accepted_hits_plus.bw"
MINUS=$ABRV_NAME/$ABRV_NAME"_accepted_hits_minus.bw"
FINAL=$ABRV_NAME/$ABRV_NAME"_accepted_hits_sorted.bw"

#if [ $REFERENCE = "NCBI" ]; then #add "Chr" to bam file
# sort the bam file and rename variable
samtools sort $CUFFIN $SBAM
SBAM="$SBAM"".bam"
#switch bam to bed to facilitate the plus/minus strand split
bamToBed -i $SBAM -split > accepted_hits.bed
awk '{if($6=="+") print}' accepted_hits.bed | sort -k1,1 | bedItemOverlapCount mm10 -chromSize=$CHROMSIZE stdin | sort -k1,1 -k2,2n > accepted_hits.plus.bedGraph
awk '{if($6=="-") print}' accepted_hits.bed | sort -k1,1 | bedItemOverlapCount mm10 -chromSize=$CHROMSIZE stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > accepted_hits.minus.bedGraph
bedGraphToBigWig accepted_hits.plus.bedGraph $CHROMSIZE $PLUS
bedGraphToBigWig accepted_hits.minus.bedGraph $CHROMSIZE $MINUS

genomeCoverageBed -bg -ibam $SBAM -g $CHROMSIZE > sorted_hits.bedgraph
bedClip -verbose=2 sorted_hits.bedgraph $CHROMSIZE clipped.bedgraph

bedGraphToBigWig clipped.bedgraph $CHROMSIZE $FINAL

rm accepted_hits.bed accepted_hits.plus.bedGraph accepted_hits.minus.bedGraph
if [ -d BIGWIGS ]; then echo ""; else mkdir BIGWIGS; fi
mv $FINAL BIGWIGS/
mv $MINUS BIGWIGS/
mv $PLUS BIGWIGS/
fi



echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`

#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID