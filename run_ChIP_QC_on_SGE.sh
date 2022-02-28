#!/bin/bash
#
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 1

#########################################################################
####################### DEFAULTS GO HERE ################################
#########################################################################

COMPRESS="y"; DATE=`date`; NAME="J. Gordon"; PATH_TO_FASTQ="./fastq"; ANALYSIS="full"

#########################################################################
####################### HELP/ USAGE      ################################
#########################################################################

DIV1=`eval printf '=%.0s' {1..100}`
DIV2=`eval printf '=%.0s' {1..25}`
usage ()
{
printf "%s\n" "" "$DIV1"
printf "%s\t" "$DIV2" "run_ChIP_QC_on_SGE.sh" "" "" "" "$DIV2"
printf "%s\n" "" "$DIV1" ""
exit 1;
}

#########################################################################
####################### Flag / Options      #############################
#########################################################################
while getopts "A:F:G:u:h:" opt; do
    case $opt in
        A) ANALYSIS=$OPTARG;;
        F) FASTQ=$OPTARG;;
        G) GENOME=$OPTARG;;
        u) NAME=$OPTARG;;
        h) usage;;
        :) usage;;
        *) usage;;
    esac
done
shift $(($OPTIND -1))
#if [[ $ANALYSIS eq "fastqc" ]]; then run fastqc; elif [[ $ANALYSIS eq "bamcorr" ]] then run bamCorrelate; elif [[ $ANALYSIS eq "IDR" ]]; then run IDR;

##########################################################
###############  FASTQC  #################################
##########################################################
if [[ -d $PATH_TO_FASTQ ]]; then FASTQS=`ls $PATH_TO_FASTQ`; echo $FASTQS; fi
mkdir fastqc
# fastqc
# md5sum for GEO before archiving


##########################################################
###############  BAMCorrelate   ##########################
##########################################################

# options for species
readarray BAMFILES < list_of_bams.txt
LABELS=
BAMS=${BAMFILES[@]}

bamCorrelate BED-file --BED RefSeq_Genes.bed --bamfiles $BAMS --labels --binSize 1000 --corMethod spearman -f 200 --colorMap Blues --zMin 0.5 --zMax 1 -o correlation_spearman.pdf

bamCorrelate BED-file --BED ../hg38_GENCODE_v22_GENES_plus_50kb_UCSC.bed --bamfiles MCF10A_ctrl_H3K27AC_nosizesel_ACAGTG_L008_R1_cat000_sort_bowtie2_hg38.bam  MCF10A_H3K27AC_R1_06032013_sorted.bam MCF10A_H3K27AC_R2_03112014_sorted.bam MCF10A_ctrl_input_nosizesel_CGATGT_L008_R1_cat000_sort_bowtie2_hg38.bam  MCF10A_input_R1_03182013_sorted.bam MCF10A_input_R1_06202013_sorted.bam --labels H3K27AC_nosize H3K27AC_REP1 H3K27AC_REP2 Nosize_Input Input_03182013 Input_06202013 --binSize 1000 --corMethod spearman -f 200 --colorMap Blues --zMin 0.5 --zMax 1 -o correlation_spearman.pdf

##########################################################
###############  IDR            ##########################
##########################################################


##########################################################
###############  bamFingerprint/ FriP ####################
##########################################################


##########################################################
###############  Dispersion from Peak summits ############
##########################################################

# calculate average distance of peak to motif.bed if provided


