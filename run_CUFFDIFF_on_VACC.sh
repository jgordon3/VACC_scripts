#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=30:00:00
#PBS -N CUFFDIFF
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

WD=/users/c/t/ctye/scratch/Runx2_treated/TOPHAT_using_no_novel_junctions/


# REQUIRES: WD(with TOPHAT and CUFFLINKS output sould have 2 folders one named "CELL_CONDITION_REP" and ""CELL_CONDITION_REP_cuffout")
# REQUIRES: SPECIES ("mm10" or 'hg19" DEFAULTS to "mm10")
# OPTION: REFERENCE (deafults to "NCBI", FNC (is file in "CELL_CONDITION_REP" convention (Y/N)), LABELS
# Don't change the following unless absolutely needed
if [ -d "$WD" ]; then cd $WD; else echo "Could not find working directory with TOPHAT and CUFFLNK files"; exit 1; fi
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
else echo "REFERENCE: $REFERENCE COULD NOT BE SET"
fi

echo "GENERATE_ASSEMBLIES_FOR_CUFFMERGE-----------"

#FILES TO PROCESS



FILES="./*_cuffout"
for f in $FILES; do
    F2=$(echo $f | awk -F '/' '{print $NF}')
#USE FOLDER NAMES FOR LABEL NAMES
    F3=$(echo $F2 | grep -E '*_cuffout'); echo "$F3" >> names.txt
    CELLS=$(echo $F3 | awk -F '_' '{print $1}')
    CONDS=$(echo $F3 | awk -F '_' '{print $2}')
    REPS=$(echo $F3 | awk -F '_' '{print $3}')
    echo "$CONDS" >> conds.txt
    echo "$REPS" >> reps.txt


echo "./$F2/transcripts.gtf" >> assemblies.txt
done
echo "GENERATED ASSEMBLIES FOR CUFFMERGE:"

cat assemblies.txt
NAMESCUFFOUNT=$(cat names.txt)
NAMES=$(cat minuscuffoutnames.txt)
UCONDS=$(cat conds.txt | sort | uniq |sed ':a;N;$!ba;s/\n/ /g')
UCONDSC=$(cat conds.txt | sort | uniq |sed ':a;N;$!ba;s/\n/,/g')
UREPS=$(cat reps.txt| sort | uniq | sed ':a;N;$!ba;s/\n/ /g')
# sed explanation create a label via :a
# append the current and next line to the pattern space via N
#if we are before the last line, branch to the created label $!ba ($! means not to do it on the last line (as there should be one final newline)).
#finally the substitution replaces every newline with a space on the pattern space (which is the whole file).

echo "GROUPS ARE: ""$UCONDS"
echo "REPS ARE: ""$UREPS"

echo "------------------------CUFFMERGE-----------------------------------"
cuffmerge -g $GENES -s $GENOMEFA assemblies.txt

echo "------------------------CUFFDIFF-----------------------------------"

CDOUT="cuffdiff_reps_out"
CDLABELS="$UCONDSC"
CDMERGE="merged_asm/merged.gtf"


#COUNT NUMBER OF UNIQUE CONDITIONS AND REPLICATES
#VAR=( $UCONDS ); COUNTCONDS=$(echo ${#VAR[@]})
#VAR=( $UREPS ); COUNTREPS=$(echo ${#VAR[@]})

for f in $UCONDS; do
    for f1 in $UREPS; do
    printf "./$CELLS" >> labels.txt; printf "_$f" >> labels.txt; printf "_$f1" >> labels.txt; printf "/accepted_hits.bam," >> labels.txt
    done
printf "\n" >> labels.txt
done
LIST=$(sed 's/,$/ /' labels.txt)
rm names.txt conds.txt reps.txt labels.txt
echo "Labels are $LIST"

echo "cuffdiff -o $CDOUT -b $GENOMEFA -L $CDLABELS -u $CDMERGE $LIST"
cuffdiff -o $CDOUT -b $GENOMEFA -L $CDLABELS -u $CDMERGE $LIST

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`

#UCONDS=$(cat conds.txt | sort | uniq |sed ':a;N;$!ba;s/\n/,/g')
#UREPS=$(cat reps.txt| sort | uniq | sed ':a;N;$!ba;s/\n/,/g')

#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID