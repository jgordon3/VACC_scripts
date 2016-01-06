#!/bin/bash
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 4

#########################################################################
####################### DEFAULTS GO HERE ################################
#########################################################################

BOWTIE_VERSION=`bowtie2 --version | head -n1`; COMPRESS="y"; DATE=`date`; FILTER="y"; GENOME="hg38";
NAME="J. Gordon"; RECIPE="SE"

#########################################################################
####################### HELP/ USAGE      ################################
#########################################################################

DIV1=`eval printf '=%.0s' {1..100}`
DIV2=`eval printf '=%.0s' {1..25}`
usage ()
{
printf "%s\n" "" "$DIV1"
printf "%s\t" "$DIV2" "run_Bowtie2_on_SGE.sh" "" "" "" "$DIV2"
printf "%s\n" "" "$DIV1" ""
printf "%s\n" "run_Bowtie2_on_SGE.sh" ""
printf "%s\n" "REQUIRES (NON OPTIONAL): A fastq file to align."
printf "%s\n" "This can be providied by the -F flag with a path/to/a/fastq"
printf "%s\n" "REQUIRES: A reference genome"
printf "%s\n" "The script will default to hg38 or the -G flag can be used. Options are hg38, hg19, mm9, mm10" ""
printf "%s\n" "OPTIONS:" "-t: Trim fastq using Trimmomatic to cut Illumina TruSeq adapters and filter for quality"
printf "%s\n" "OPTIONS:" "This option generates an extra (timmed) fastq file"
printf "%s\n" "OPTIONS:" "-c: Compress outputs fastq(s) to gz then tar"
exit 1;
}

#########################################################################
####################### Flag / Options      #############################
#########################################################################

while getopts "F:c:G:t:u:q:R:h:" opt; do
    case $opt in
        F) FASTQ=$OPTARG;;
        c) COMPRESS=$OPTARG;;
        G) GENOME=$OPTARG;;
        t) FILTER=$OPTARG;;
        u) NAME=$OPTARG;;
        q) FASTQC=$OPTARG;;
        R) RECIPE=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

###################### CHECK FASTQ(S) #####################################
if [[ -e $FASTQ ]]; then FASTQ=$FASTQ; else echo "Need valid path to fastq or fastq.gz file"; usage; fi
#if [[ $RECIPE == "PE" ]]; then check for pair ; fi

#CHECK FOR GZIP
EXT=${FASTQ##*.}
if [ $EXT = "gz" ]; then gunzip $FASTQ; FASTQ=${FASTQ%.gz}; fi

# FILE NAME
ABRV_NAME=${FASTQ%.fastq}



# GENOME ASSIGNMENT
if [ $GENOME == "hg38" ]; then INDEX=/slipstream/galaxy/data/hg38/hg38canon/bowtie2_index/hg38canon
        elif [ $GENOME == "hg19" ]; then INDEX=/slipstream/galaxy/data/hg19/hg19canon/bowtie2_index
        elif [ $GENOME == "mm10" ]; then INDEX=/slipstream/galaxy/data/mm10/bowtie2_index_canon
        elif [ $GENOME == "mm9" ]; then INDEX=/slipstream/galaxy/data/mm9/bowtie2_index_canon
else echo "could not find specified index"; usage;
fi

# OPTIONAL PREFILTERING USING TRIMMOMATIC
TRIMMOMATIC_DIR=/slipstream/home/jonathan/bin/Trimmomatic-0.35
TRIMMOMATIC_ADAPTERS_DIR=/slipstream/home/jonathan/bin/Trimmomatic-0.35/adapters
# TO DO: give adapter file options - Paired end and single end
if [[ $FILTER = [Yy] ]]; then
    printf "%s\n" "$DIVIDER" "" "$FASTQ was filtered with Trimmomatic using the following parameters:" ""
    FILT_FASTQ=$ABRV_NAME"_trim_filt.fastq";
    java -Djava.io.tmpdir=/slipstream/home/jonathan/bin/tmp -jar $TRIMMOMATIC_DIR/trimmomatic-0.35.jar SE -threads 4 -phred33 $FASTQ $FILT_FASTQ ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS_DIR/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36;
    else  printf "%s\n" "" "No prefiltering was selected" ""; FILT_FASTQ=$FASTQ
fi

######################################################################################
#################### run BOWTIE2 here ################################################
######################################################################################

printf "%s\n" "$DIVIDER" ""  "RUNNING BOWTIE2 ON: $FILT_FASTQ" ""
OUTPUT_NAME=$ABRV_NAME"_BAM.bowtie"
bowtie2 -p 4 -x $INDEX -U $FILT_FASTQ -S $OUTPUT_NAME


######################################################################################
#################### SAMTools Sort and Index #########################################
######################################################################################

printf "%s\n" "$DIVIDER" "Sorting the Bowtie output and indexing" ""
SORTED_NAME=$ABRV_NAME"_sort_bowtie2_"$GENOME
samtools view -Su $OUTPUT_NAME | samtools sort - $SORTED_NAME
SORTED_NAME=$SORTED_NAME".bam"
samtools index $SORTED_NAME
SAM_INDEX=$SORTED_NAME".bai"
# remove unsorted bowtie file
rm $OUTPUT_NAME

# Compress fastqs and sorted sam files

if [[ $COMPRESS = [Yy] ]]; then
    printf "%s\n" "Compressing the output(s)......" ""
    # Compress fastqs (2)
    gzip $FILT_FASTQ
    FASTQ_TAR=$FILT_FASTQ
    if [[ "$FASTQ" ! = "$FILT_FASTQ" ]]; then # this is true if FILTER=y
            gzip $FASTQ;
            FASTQ=$FASTQ".gz"; FILT_FASTQ=$FILT_FASTQ".gz"; FASTQ_TAR=$ABRV_NAME"_fastqs.tar"
            tar -cvf $FASTQ_TAR $FASTQ $FILT_FASTQ --remove-files
    fi

    # Compress bam files
    gzip $SORTED_NAME $SAM_INDEX
    SORTED_NAME=$SORTED_NAME".gz"; SAM_INDEX=$SAM_INDEX".gz"; SAM_TAR=$ABRV_NAME"_bamfiles.tar"
    tar -cvf $SAM_TAR $SORTED_NAME $SAM_INDEX --remove-files
    printf "%s\n" "GENERATED TWO ARCHIVES: $FASTQ_TAR and $SAM_TAR" "" "$DIVIDER"
    else printf "%s\n" "" "OUTPUTS WERE NOT COMPRESED" "" "$DIVIDER"
fi

# HEADER
printf "%s\n" "$DIVIDER" "" "$ABRV_NAME" "" "ALIGNED TO $GENOME" "" "BOWTIE VERSION: $BOWTIE_VERSION" "" "USER: $NAME" "" "DATE: $DATE" ""

# Rename standard out file
README=$ABRV_NAME"_Bowtie2_align_README.txt"
mv *.o$JOB_ID $README
rm *.po$JOB_ID




#PASSED=$JOB_ID
#qsub -v PASSED="$JOB_ID" ~/test.sh