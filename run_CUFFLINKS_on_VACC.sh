#!/bin/bash
#PBS -l nodes=1:ppn=8,pmem=1.5gb,pvmem=2gb
#PBS -l walltime=08:00:00
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea













#REQUIRES: SPECIES,INPUT (TOPHAT folder containing "accepted_hits")
#OPTIONS: WD (WD will be set to INPUT foleder by default), REFERENCE

#WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE
if [ -z $WD ]; then WD=$INPUT; cd $WD; else echo "WORKING DIRECTORY WAS SET TO $WD"; cd $WD; fi
if [ -d $INPUT ]; then echo "INPUT DIRECTORY IS $INPUT"; else echo "NO INPUT DIRECTORY: QUIT"; exit 1; fi
rm $WAIT

#---------LOCATION_OF_GENE_REFERENCE_FILES---------
if [ -z $SPECIES ]; then SPECIES="mm10"; echo "GENOME HAS BEEN SET TO $SPECIES BY DEFAULT"; else printf "."; fi
if [ -z $REFERENCE ]; then REFERENCE="UCSC"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else printf "."; fi
if [ $SPECIES = "mm10" ]; then FULLSPECIES="Mus_musculus"; elif [ $SPECIES = "hg19" ]; then FULLSPECIES="Homo_sapiens"; else echo "REFERENCE: $REFERENCE COULD NOT BE SET"; fi

INDEXFOLDER=/users/j/a/jargordo/scratch/indexes/$FULLSPECIES

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
#----------FILES_TO_PROCESS-----------------------

CUFFIN=$INPUT"/accepted_hits.bam"
FNAME=$(echo $INPUT | awk -F '/' '{print toupper($NF)}' | awk -F "." '{print $1}')
CUFFOUT="$FNAME""_cuffout"

# BUILD CUFFLINKS STATEMENT
QSUB="cufflinks -g $GENES"

if [ -z $LIBTYPE ]; then echo ""; else QSUB="$QSUB --library-type=$LIBTYPE "; fi

QSUB="$QSUB -o $CUFFOUT $CUFFIN"  

echo "------------------------CUFFLINKS-----------------------------------"
echo `$QSUB`

#CREATE assemblies.txt
#echo "$CUFFOUT""/transcripts.gtf" >> assemblies.txt

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`



#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID