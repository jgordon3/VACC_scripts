#!/bin/bash
#$ -N STAR
#$ -cwd
#$ -j y
# Request [num] amount of [type] nodes
#$ -pe threads 10
# STAR path: /slipstream/galaxy/production/galaxy-dist/tools/star/STAR

######################## DEFAULTS #############################################

DATE=`date`; GENOME="hg38"; NAME="J. Gordon"; RECIPE="SE"; MODE="RNA";

######################## USAGE ###############################################

usage ()
{
printf "%s\n" "" "$DIV1"
printf "%s\t" "$DIV2" "run_STAR_on_SGE.sh" "" "" "" "$DIV2"
printf "%s\n" "" "$DIV1" ""
printf "%s\n" "run_STAR_on_SGE.sh" ""
printf "%s\n" "REQUIRES (NON OPTIONAL): A fastq file to align."
printf "%s\n" "This can be providied by the -F flag with a path/to/a/fastq"
printf "%s\n" "REQUIRES: A reference genome"
printf "%s\n" "OPTIONS: -G: GENOME can be hg38 (default), hg19, mm9, mm10"
printf "%s\n" "OPTIONS: -R: RECIPE can be single end SE (default) or paired end PE"
printf "%s\n" "OPTIONS: -M: MODE can be RNA (defaut) or CHIP"
exit 1;
}

######################## FLAGS ##################################################

while getopts "F:G:t:u:R:M:h:" opt; do
    case $opt in
        F) FASTQ=$OPTARG;;
        G) GENOME=$OPTARG;;
        t) FILTER=$OPTARG;;
        u) NAME=$OPTARG;;
        R) RECIPE=$OPTARG;;
        M) MODE=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))

###################### CHECKS #####################################
# FASTQ
FASTQ=$(readlink -f $FASTQ)
if [ ! -f "$FASTQ" ]; then echo "Need valid path to fastq or fastq.gz file"; fi

# PAIRED
if [[ $RECIPE == "PE" ]]; then FASTQ2=${FASTQ/_R1_/_R2_};
    if [ ! -f $FASTQ2 ]; then echo "Paired end flag used (-R PE) and a matching pair could not be found."; echo 'Check that file name is in: "Some_fastq_R2_000.fastq" format'; usage;
    fi
fi

# GZIP
FASTQEXT=${FASTQ##*.}
if [ $FASTQEXT = "gz" ]; then gunzip $FASTQ; FASTQ=${FASTQ%.gz}; fi
FASTQ2EXT=${FASTQ2##*.}
if [ -f $FASTQ2 ] && [[ $FASTQ2EXT = "gz" ]]; then gunzip $FASTQ2; FASTQ2=${FASTQ2%.gz}; fi

# OUTPUT FOLDER
if [ ! -d "./BAMS" ]; then mkdir BAMS; fi

###################### FILE NAMES #####################################
BASE=`basename $FASTQ`
ABRV_NAME=${BASE%.fastq}
mkdir $ABRV_NAME
cd $ABRV_NAME

# GENOME ASSIGNMENT
if [ $GENOME == "hg38" ]; then INDEX=/slipstream/home/jonathan/MALIK_NATURE_COMMUN_2019/RNAseq/hg38_index
elif [ $GENOME == "hg19" ]; then INDEX=/slipstream/galaxy/data/hg19/star_index
elif [ $GENOME == "mm10" ]; then INDEX=/slipstream/galaxy/data/mm10/star_index
elif [ $GENOME == "mm9" ]; then INDEX=/slipstream/galaxy/data/mm9/star_index
else echo "could not find specified index"; usage;
fi
$GTF=/slipstream/home/jonathan/MALIK_NATURE_COMMUN_2019/RNAseq/gencode.v38.annotation.gtf

##################### RUN STAR ##########################################################

STAR \
--runThreadN 10 \
--genomeDir $INDEX \
--sjdbGTFfile $GTF \
--runThreadN 12 \
--readFilesIn $FASTQ $FASTQ2 \
--readFilesCommand gunzip -c \
--runDirPerm User_RWX \
--outFilterType BySJout\
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--outFileNamePrefix $ABRV_NAME \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--twopassMode Basic





# TROUBLESHOOT ######
#echo "FASTQ: $FASTQ"; echo "FASTQ2: $FASTQ2"; echo "FASTQEXT: $FASTQEXT"; echo "FASTQ2EXT: $FASTQ2EXT";
#echo "BASE: $BASE"; echo "ABRV_NAME: $ABRV_NAME"; echo "OUTPUT_NAME: $OUTPUT_NAME";
#echo "RECIPE: $RECIPE"; echo "GENOME: $GENOME"; echo "MODE: $MODE";




