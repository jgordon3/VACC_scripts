#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=24:00:00
#PBS -N CT_TOPHAT
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

#REQUIRES: SPECIES
#OPTIONS:WD (working directory), REFERENCE (will be set to NCBI by default), LIBTYPE (Library type options for TopHat and Cufflinks: default is unstranded)
# SEE: Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., … Pachter, L. (2012).
# Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.
# Nature Protocols, 7(3), 562–78. doi:10.1038/nprot.2012.016
# don't change this script unless absolutely needed


# SETTING DIRECTORIES
DATE=`date +%m_%d_%Y`
if [ -z $WD ]; then WD=$HOME/scratch/"TOPHAT_"$DATE; else echo "working directory is $WD"; fi
if [ -d "$WD" ]; then cd $WD; else mkdir $WD; cd $WD; fi

#---------LOCATION_OF_GENE_REFERENCE_FILES---------
#need to finish
if [ -z $SPECIES ]; then SPECIES="mm10"; echo "GENOME HAS BEEN SET TO $SPECIES BY DEFAULT"; else echo "GENOME WAS SET TO $SPECIES"; fi
if [ -z $REFERENCE ]; then REFERENCE="UCSC"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else echo "REFERENCE WAS SET TO $REFERENCE"; fi
if [ $SPECIES = "mm10" ]; then FULLSPECIES="Mus_musculus"; BUILD="GRCm38";
    elif [ $SPECIES = "hg19" ]; then FULLSPECIES="Homo_sapiens"; BUILD="build37.2";
else echo "REFERENCE: $REFERENCE COULD NOT BE SET"; exit 1;
fi
#REFERENCE
if [ $REFERENCE = "UCSC" ]; then
    echo "ALIGNMENT TO: IGENOMES UCSC $FULLSPECIES $SPECIES"
    GENES=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/UCSC/$SPECIES/Annotation/Genes/genes.gtf
    GENOME=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/UCSC/$SPECIES/Sequence/Bowtie2Index/genome
    GENOMEFA=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/UCSC/$SPECIES/Sequence/Bowtie2Index/genome.fa
    CHROMSIZE=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/UCSC/$SPECIES/Annotation/Genes/ChromInfo.txt
elif [ $REFERENCE = "NCBI" ]; then
    echo "ALIGNMENT TO: IGENOMES $FULLSPECIES  $SPECIES GRCm38 (NCBI)"
    GENES=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/NCBI/$BUILD/Annotation/Genes/genes.gtf
    GENOME=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/NCBI/$BUILD/Sequence/Bowtie2Index/genome
    GENOMEFA=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/NCBI/$BUILD/Sequence/WholeGenomeFasta/genome.fa
    CHROMSIZE=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES/NCBI/$BUILD/GenomeStudio/$FULLSPECIES/NCBI-$BUILD/ChromInfo.txt
else echo "REFERENCE: $REFERENCE COULD NOT BE SET"
fi

# DEFINE INPUTS/OUTPUTS
# CHECK FOR GZIP
GZIP=$(echo $PATHTOFASTQ | awk -F '.' '{print $NF}')
if [ "$GZIP" = "gz" ]; then gunzip $PATHTOFASTQ; PATHTOFASTQ=$(echo $PATHTOFASTQ| rev | cut -d. -f2- | rev); fi

# CHECK FOR LIBTYPE
if [ -z $LIBTYPE ]; then LIBTYPE="fr-firststrand"; fi


#FILE NAMING STUFF
FILE=$(echo $PATHTOFASTQ | awk -F '/' '{print $NF}')
ABRV_NAME=$(echo "$FILE" | awk -F '.' '{print $1}')

CUFFOUT="$WD/$ABRV_NAME/$ABRV_NAME""_cuffout"

echo "INPUT FILE: $FILE"
echo "OUTPUT FILE(S): $CUFFOUT"

# SET ANCHOR LENGTH
if [ -z $ALENGTH ]; then ALENGTH=5; fi


#STUFF ACTUALLY GETS DONE HERE
echo "TOPHAT CALL: tophat -G $GENES -o $ABRV_NAME -a $ALENGTH --library-type=$LIBTYPE $GENOME $PATHTOFASTQ"

echo "START_OF_TOPHAT_OUTPUT-----------------------------------"
tophat -G $GENES -o $ABRV_NAME -a $ALENGTH --library-type=$LIBTYPE $GENOME $PATHTOFASTQ



#SORT --> BAM --> BEGRAPH --> BIGWIG
CUFFIN=$WD/$ABRV_NAME/"accepted_hits.bam"
SBAM=$WD/$ABRV_NAME/"accepted_hits_sorted"
PLUS=$WD/$ABRV_NAME/$ABRV_NAME"_accepted_hits_plus.bw"
MINUS=$WD/$ABRV_NAME/$ABRV_NAME"_accepted_hits_minus.bw"
FINAL=$WD/$ABRV_NAME/$ABRV_NAME"_accepted_hits_sorted.bw"

#if [ $REFERENCE = "NCBI" ]; then #add "Chr" to bam file
# sort the bam file and rename variable
samtools sort $CUFFIN $SBAM
SBAM="$SBAM"".bam"
#switch bam to bed to facilitate the plus/minus strand split
bamToBed -i $SBAM -split > accepted_hits.bed
awk '{if($6=="+") print}' accepted_hits.bed | sort -k1,1 | bedItemOverlapCount $SPECIES -chromSize=$CHROMSIZE stdin | sort -k1,1 -k2,2n > accepted_hits.plus.bedGraph
awk '{if($6=="-") print}' accepted_hits.bed | sort -k1,1 | bedItemOverlapCount $SPECIES -chromSize=$CHROMSIZE stdin | sort -k1,1 -k2,2n | awk '{OFS="\t"; print $1,$2,$3,"-"$4}' > accepted_hits.minus.bedGraph
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


echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`

#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID