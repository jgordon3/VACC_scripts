#!/bin/bash
#$ -N HISAT2_TOPHAT2
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 4

#GENERAL STUFF
DIVIDER=`eval printf '=%.0s' {1..100}`
TOPHAT_VERSION=`tophat2 --version`
DATE=`date`


# SEE: Trapnell, C., Roberts, A., Goff, L., Pertea, G., Kim, D., Kelley, D. R., … Pachter, L. (2012).
# Differential gene and transcript expression analysis of RNA-seq experiments with TopHat and Cufflinks.
# Nature Protocols, 7(3), 562–78. doi:10.1038/nprot.2012.016
# don't change this script unless absolutely needed

usage ()
{
case $ERROR in
    1) MESSAGE="You did not supply a fastq with the -F flag or it does not exist";;
    2) MESSAGE="Paired end flagged (-p) but no mate file";;
    3) MESSAGE="Genome was not set ... which is weird because it is supposed to default to hg38. Please set it with -G";;
    *) MESSAGE="Help file activated or Random Error: resubmit";;
esac
printf "%s\n" "$DIVIDER" "" "$MESSAGE" "" "This is script is now submitted via qsub with the command" "";
printf "%s\n" "qsub run_tophat2_hisat2_on_SGE.sh -F /path/to/fastq" ""
printf "%s\n" "****  optional flags: ******** " "" "-p: Paired end data and will search for matched read (R2) unless a path to the mate is supplied. Default assumes no pair" "-G: Genome for alignment. Default is hg38. Some other options are hg19, mm9, mm10 ... etc" "-t: [y/n]: filter fastq with Trimmomatic for Illumina Adapters, quality, and length. Default is y" "-u: User info for README file" "-q: [y/n] Submit fastqc script for raw and/or trimmed fastq files. Default is n." "$DIVIDER" ""
exit 1;
}

#set flags here

while getopts "F:pG:th" opt; do
        case $opt in
        F) FASTQ=$OPTARG;;          #fastq is required
        p) PAIREDEND=1;;
        G) GENOME=$OPTARG;;
        t) FILTER="y";;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
        esac
    done
shift $(($OPTIND -1))

printf "%s\n" "$DIVIDER" ""

if [ -e $FASTQ ]; then printf "%s\n" "Fastq: $FASTQ"; else ERROR=1; usage; fi
if [[ $HELP = 1 ]]; then usage; fi
if [ -z $GENOME ]; then GENOME="mm10"; printf "%s\n" "Genome set by default to: $GENOME"; else GENOME=$GENOME; printf "%s\n" "Genome set by flag (-G) to: $GENOME"; fi

#CHECK FOR PAIRED FILE
if [[ $PAIREDEND = 1 ]]; then FASTQ2=${FASTQ//_R1_/_R2_}; TYPE="PE100"; else TYPE="SE100"; fi  #fix this to be more general # check to see if there is more than one R1 in file name
if [[ $PAIREDEND = 1 ]] && [ ! -e $FASTQ2 ]; then ERROR=2; usage; else printf "%s\n" "Paired fastq was found: $FASTQ2"; fi


# CHECK FOR GZIP
EXT=${FASTQ##*.}
if [ "$EXT" = "gz" ]; then gunzip $FASTQ; gunzip $FASTQ2; FASTQ=${FASTQ%.gz}; FASTQ2=${FASTQ2%.gz}; printf "%s\n" "Fastq is gzipped ..... unzipping"; fi
printf "%s\n" "$DIVIDER" ""

ABRV_NAME=${FASTQ%.fastq}
OUTPUT_NAME=$ABRV_NAME"_TOPHAT_OUT"; mkdir $OUTPUT_NAME;

# OPTIONAL PREFILTERING USING TRIMMOMATIC
TRIMMOMATIC_DIR=/slipstream/home/jonathan/Trimmomatic-0.33
TRIMMOMATIC_ADAPTERS_DIR=/slipstream/home/jonathan/Trimmomatic-0.33/adapters
    # TO DO: give adapter file options
if [[ $FILTER = [Yy] ]]; then
        printf "%s\n" "$DIVIDER" "" "$FASTQ was filtered with Trimmomatic using the following parameters:" ""
        FILT_FASTQ=$ABRV_NAME"_trim_filt.fastq";
        java -jar $TRIMMOMATIC_DIR/trimmomatic-0.33.jar SE -threads 4 -phred33 $FASTQ $FILT_FASTQ ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS_DIR/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36;
        else  printf "%s\n" "" "No prefiltering was selected" ""; FILT_FASTQ=$FASTQ
fi

# SET GENOME INDEX AND GTF #find GTF paths
if [ "$GENOME" = "hg38" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/hg38/hg38canon/bowtie2_index/hg38canon;
    elif [ "$GENOME" = "hg19" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/hg19canon/bowtie2_index_canon/hg19canon;
    elif [ "$GENOME" = "mm10" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/mm10/bowtie2_index_canon/mm10canon;
    elif [ "$GENOME" = "mm9" ]; then BOWTIE_INDEX_PATH=/slipstream/galaxy/data/mm9/bowtie2_index_canon/mm9canon;
    else ERROR=3; usage;
fi

#CONDITIONAL OR CUSTOM GTF
REF_GTF=REFGENE_MM10_REFFLAT_UCSC.GTF

tophat2 -p 4 --output-dir $OUTPUT_NAME -G $REF_GTF $BOWTIE_INDEX_PATH $FASTQ $FASTQ2

#RENAME FILES: accepted_hits.bam  align_summary.txt  deletions.bed  insertions.bed  junctions.bed  logs  prep_reads.info  unmapped.bam
for f in $OUTPUT_NAME/*; do FILE=`basename $f`; N=$ABRV_NAME"_"$FILE;  mv $f $OUTPUT_NAME/$N; done

BAM_FILE=$ABRV_NAME"_accepted_hits.bam"

# CREATE SUMMARY with alignments and other parameters

SUMMARY=$ABRV_NAME"_align_summary.txt"
printf "%s\n" "" "FASTQ: $FASTQ" "" "GENOME: $GENOME" "" "GTF: REF_GTF" "" "TYPE: $TYPE" "" "" "VERSION: $TOPHAT_VERSION" >> $SUMMARY

#REMOVE STDOUT - Launch a clean script?

#LAUNCH HTSEQ-COUNT OR CUFFDIFF (flag this?)
qsub ~/run_htseq_on_SGE.sh -D $OUTPUT_NAME -G $REF_GTF

#LAUNCH BIGWIG CONVERSION
