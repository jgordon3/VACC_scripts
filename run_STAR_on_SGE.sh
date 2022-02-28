#!/bin/bash
#$ -N STAR
#$ -cwd
#$ -j y
# Request [num] amount of [type] nodes
#$ -pe threads 8
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
OUTPUT_NAME=$ABRV_NAME"_STAR_ALIGNMENT.bam"

# GENOME ASSIGNMENT
if [ $GENOME == "hg38" ]; then INDEX=/slipstream/galaxy/data/hg38/star_index
elif [ $GENOME == "hg19" ]; then INDEX=/slipstream/galaxy/data/hg19/star_index
elif [ $GENOME == "mm10" ]; then INDEX=/slipstream/galaxy/data/mm10/star_index
elif [ $GENOME == "mm9" ]; then INDEX=/slipstream/galaxy/data/mm9/star_index
else echo "could not find specified index"; usage;
fi


##################### RUN STAR ##########################################################

if [ $MODE == "RNA" ]; then
STAR --genomeLoad NoSharedMemory --genomeDir $INDEX --readFilesIn $FASTQ $FASTQ2 --runThreadN 8 --seedSearchStartLmax 30 --outTmpDir $ABRV_NAME --outFileNamePrefix $ABRV_NAME --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > $OUTPUT_NAME
elif [ $MODE == "CHIP" ]; then
STAR --genomeLoad NoSharedMemory --genomeDir $INDEX --readFilesIn $FASTQ $FASTQ2 --runThreadN 8 --alignIntronMax 1 --seedSearchStartLmax 30 --outTmpDir $ABRV_NAME --outFileNamePrefix $ABRV_NAME --outSAMtype BAM SortedByCoordinate --outStd BAM_SortedByCoordinate > $OUTPUT_NAME
else usage
fi

##################### CLEAN UP ##########################################################
# GZIP
if [ $FASTQEXT = "gz" ]; then gzip $FASTQ; fi
if [ -f $FASTQ2 ] && [[ $FASTQ2EXT = "gz" ]]; then gzip $FASTQ2; fi

# REMOVE SORT BAMS
rm -r $ABRV_NAME

# MOVE BAM
mv $OUTPUT_NAME ./BAMS/

# REMOVE UNNECCESSARY LOGS
rm $ABRV_NAME"Log.out" $ABRV_NAME"Log.progress.out" $ABRV_NAME"Log.std.out" $ABRV_NAME"SJ.out.tab"

# CREATE SUMMARY
ALIGN_SUM="STAR_alignment_summaries.txt"
COMB_LOG="STAR_combined_logs.txt"
STAR_LOG="$ABRV_NAME""Log.final.out"
if [ ! -f $ALIGN_SUM ]; then touch $ALIGN_SUM;
printf "%s\t" "File" >> $ALIGN_SUM;
while IFS='' read -r line || [[ -n "$line" ]]; do
new=$(echo $line | awk '{$1=$1}1');
split=$(echo $new | awk -F\| '{print $1}');
printf "%s\t" ${split// /_} >> $ALIGN_SUM;
done < "$STAR_LOG"
printf "\n" >> $ALIGN_SUM
fi

printf "%s\t" "$ABRV_NAME" >> $ALIGN_SUM;
while IFS='' read -r line || [[ -n "$line" ]]; do
new=$(echo $line | awk '{$1=$1}1');
split=$(echo $new | awk -F\| '{print $2}');
printf "%s\t" ${split// /} >> $ALIGN_SUM;
done < "$STAR_LOG"
printf "\n" >> $ALIGN_SUM

echo "$ABRV_NAME" >> $COMB_LOG
cat $STAR_LOG >> $COMB_LOG
printf "\n"
rm $STAR_LOG

# TROUBLESHOOT ######
#echo "FASTQ: $FASTQ"; echo "FASTQ2: $FASTQ2"; echo "FASTQEXT: $FASTQEXT"; echo "FASTQ2EXT: $FASTQ2EXT";
#echo "BASE: $BASE"; echo "ABRV_NAME: $ABRV_NAME"; echo "OUTPUT_NAME: $OUTPUT_NAME";
#echo "RECIPE: $RECIPE"; echo "GENOME: $GENOME"; echo "MODE: $MODE";




