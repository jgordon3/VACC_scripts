#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=30:00:00
#PBS -N CUFFDIFF_mBMSC
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

CUFFOUTDIR=$WD/"*_cuffout"
eval CUFFOUTDIR=$CUFFOUTDIR

# REQUIRES: WD(with TOPHAT and CUFFLINKS output sould have 2 folders one named "CELL_CONDITION_REP" and ""CELL_CONDITION_REP_cuffout")
# REQUIRES: SPECIES ("mm10" or 'hg19" DEFAULTS to "mm10")
# OPTION: REFERENCE (deafults to "NCBI", FNC (is file in "CELL_CONDITION_REP" convention (Y/N)), LABELS
# Don't change the following unless absolutely needed
if [ -d "$WD" ]; then cd $WD; else echo "Could not find working directory with TOPHAT and CUFFLNK files"; exit 1; fi
if [ -z $WD ]; then WD=$INPUT; cd $WD; else echo "WORKING DIRECTORY WAS SET TO $WD"; cd $WD; fi
#if [ -d $INPUT ]; then echo "INPUT DIRECTORY IS $INPUT"; else echo "NO INPUT DIRECTORY: QUIT"; exit 1; fi
#rm $WAIT

#---------LOCATION_OF_GENE_REFERENCE_FILES---------
if [ -z $SPECIES ]; then SPECIES="mm10"; echo "GENOME HAS BEEN SET TO $SPECIES BY DEFAULT"; else printf "."; fi
if [ -z $REFERENCE ]; then REFERENCE="UCSC"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else printf "."; fi
if [ $SPECIES = "mm10" ]; then FULLSPECIES="Mus_musculus"; elif [ $SPECIES = "hg19" ]; then FULLSPECIES="Homo_sapiens"; else echo "REFERENCE: $REFERENCE COULD NOT BE SET"; fi
if [ -z $LIBTYPE ]; then LIBTYPE="fr-unstranded"; echo "REFERENCE HAS BEEN SET TO $REFERENCE BY DEFAULT"; else printf "."; fi

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
    GENOME=$INDEXFOLDER/ENSEMBL/Sequence/Bowtie2Index/genome
    GENOMEFA=$INDEXFOLDER/ENSEMBL/Sequence/WholeGenomeFasta/genome.fa
    CHROMSIZE=$INDEXFOLDER/GenomeStudio/$FULLSPECIES/ChromInfo.txt
else echo "REFERENCE: $REFERENCE COULD NOT BE SET"
fi

echo "GENERATE_ASSEMBLIES_FOR_CUFFMERGE-----------"

#FILES TO PROCESS
rm assemblies.txt

for f in $CUFFOUTDIR; do
echo "$f/transcripts.gtf" >> assemblies.txt
done

printf "GENERATED ASSEMBLIES FOR CUFFMERGE: \n"


echo "CUFFMERGE"
cuffmerge -g $GENES -s $GENOMEFA assemblies.txt

echo "CUFFDIFF"

CDOUT="cuffdiff_out"
CDMERGE="merged_asm/merged.gtf"
CDLABELS=

echo "cuffdiff -o $CDOUT -b $GENOMEFA -L $CDLABELS -u $CDMERGE $LIST"
cuffdiff -o $CDOUT -b $GENOMEFA -L $CDLABELS -u $CDMERGE ISHIKAWA_BZA_E2_R1_R1/accepted_hits.bam,ISHIKAWA_BZA_E2_R2_R1/accepted_hits.bam ISHIKAWA_BZA_GC10_R1_R1/accepted_hits.bam,ISHIKAWA_BZA_GC10_R2_R1/accepted_hits.bam ISHIKAWA_BZA_R1_R1/accepted_hits.bam,ISHIKAWA_BZA_R2_R1/accepted_hits.bam ISHIKAWA_CONTROL_R1_R1/accepted_hits.bam,ISHIKAWA_CONTROL_R2_R1/accepted_hits.bam
echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`

#UCONDS=$(cat conds.txt | sort | uniq |sed ':a;N;$!ba;s/\n/,/g')
#UREPS=$(cat reps.txt| sort | uniq | sed ':a;N;$!ba;s/\n/,/g')

#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID