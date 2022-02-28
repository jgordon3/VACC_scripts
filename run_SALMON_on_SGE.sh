#!/bin/bash
#$ -N SALMON
#$ -cwd
#$ -j y
# Request [num] amount of [type] nodes
#$ -pe threads 8

######################## PATH #####################################
SALMONPATH="/slipstream/home/jonathan/bin/salmon/bin/"

######################## USAGE ####################################

usage ()
{
printf "%s\n" "" "$DIV1"
printf "%s\t" "$DIV2" "run_SALMON_on_SGE.sh" "" "" "" "$DIV2"
printf "%s\n" "" "$DIV1" ""
printf "%s\n" "Highly-accurate & wicked fast transcript-level quantification from RNA-seq reads using lightweight alignments"
printf "%s\n" "https://combine-lab.github.io/salmon/getting_started/"
printf "%s\n" "https://salmon.readthedocs.io/en/latest/"
printf "%s\n" "REQUIRES (NON OPTIONAL): A fastq file to align."
printf "%s\n" "can be set with the -F flag OR as the first positional argument (e.g. run_SALMON_on_SGE.sh /path/to/fastq"
printf "%s\n" "OPTIONAL FLAGS:"
printf "%s\n" "-I transcriptome index that is precompiled (default is gencode_v32_plus_novel_lncRNA_idx)"
printf "%s\n" "-R: RECIPE can be single end SE (default) or paired end PE (currently not working)"
exit 1;
}

################## OPTION FLAGS ###############################################

while getopts "F:I:t:u:R:M:h:" opt; do
    case $opt in
        F) FASTQ=$OPTARG;;
        I) INDEX=$OPTARG;;
        gc) GCBIAS=$OPTARG;;
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
#FASTQ
if [ -z $FASTQ ];then FASTQ=$1;fi
FASTQ=$(readlink -f $FASTQ)
if [ ! -f "$FASTQ" ]; then echo "Need valid path to fastq or fastq.gz file"; usage; fi
#INDEX
if [ -z $INDEX ];then INDEX="/slipstream/home/jonathan/bin/salmon/gencode_v32_plus_novel_lncRNA_idx"; fi
if [ ! -d "$INDEX" ]; then echo "Need valid path to transcriptome index"; usage; fi
# PAIRED
#if [[ $RECIPE == "PE" ]]; then FASTQ2=${FASTQ/_R1_/_R2_};
#    if [ ! -f $FASTQ2 ]; then echo "Paired end flag used (-R PE) and a matching pair could not be found."; echo 'Check that file name is in: "Some_fastq_R2_000.fastq" format'; usage;
#    fi
#fi

#GCBIAS
#if [ $GCBIAS = "n" ]; then OFLAGS=$OFLAGS; else OFLAGS=$"$OFLAGS --gcBias"; fi ..... finish this later

# GZIP
#FASTQEXT=${FASTQ##*.}
#if [ $FASTQEXT = "gz" ]; then gunzip $FASTQ; FASTQ=${FASTQ%.gz}; fi
#FASTQ2EXT=${FASTQ2##*.}
#if [ -f $FASTQ2 ] && [[ $FASTQ2EXT = "gz" ]]; then gunzip $FASTQ2; FASTQ2=${FASTQ2%.gz}; fi

# OUTPUT FOLDER
#if [ ! -d "./BAMS" ]; then mkdir BAMS; fi

###################### FILE NAMES #####################################
FASTQBASE=`basename $FASTQ`
FASTQNAME=${FASTQBASE%.fastq}
QUANTS="quants/${FASTQNAME}_quant/quant.sf"
###################### SALMON #########################################
    
salmon quant -i ${INDEX} -l A -r ${FASTQ} -p 8 --gcBias --validateMappings -o quants/${FASTQNAME}_quant

#################### MOVE QUANT FILE ####################################
while [ ! -f "$QUANTS" ]
do
  sleep 2
done
cp ${QUANTS} quants/${FASTQNAME}_quants.sf
