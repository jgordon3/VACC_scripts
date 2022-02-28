#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=04:00:00
#PBS -N MACS2
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

#READ THIS
#THINGS TO CHANGE:
# This is scipt is now submitted via the submit script but can be submited directly with the command "qsub -v BOWTIE="PATH/TO/BOWTIE/FILE",GENOME="hg19" run_MACS2_on_VACC.sh"
# nothing needs to be changed on this script
# REQUIRED: BOWTIE=path to a bowtie (SAM) file to process
# REQUIRED: GENOME=hg19 or mm10

#WORKING_DIRECTORY
DATE=`date +%m_%d_%Y`
if [ -z $WD ]; then WD=$HOME/scratch/$PBS_JOBNAME"_"$DATE; else echo "working directory is $WD"; fi
if [ -d "$WD" ]; then cd $WD; else mkdir $WD; cd $WD; fi

#SPECIES
if [ -z $GENOME ]; then GENOME="mm10"; else echo "$GENOME"; fi
if  [ $GENOME = "mm10" ]; then CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Mus_musculus/UCSC/mm10/GenomeStudio/Mus_musculus/UCSC-mm10/mm10.chrom.sizes"; GEN="mm";
elif [ $GENOME = "hg19" ]; then CHROMSIZE="/users/j/a/jargordo/scratch/indexes/Homo_sapiens/UCSC/hg19/GenomeStudio/Homo_sapiens/UCSC-hg19/hg19.chrom.sizes"; GEN="hs";
else echo "$GENOME is not set"; fi

#CHECK INPUTS AND NAME
if [ -z $BOWTIE ]; then echo "NO BOWTIE FILE!"; exit 1; else echo "TEST FILE: $BOWTIE"; fi
FILE=$(echo "$BOWTIE" | awk -F '/' '{print $NF}')
ABRV_NAME=$(echo "$FILE" | awk -F '.' '{print $1}')
NAME=$(echo "$ABRV_NAME" | awk -F '.' '{print $1 "_" $2 "_" $3}')
if [ -z $BOWTIE2 ]; then echo "NO REPLICATE FILE"; else echo "SECOND FILE FOUND: $BOWTIE2. THIS WILL COMBINE BOTH REPLICATES"; FILE2=$(echo "$BOWTIE2" | awk -F '/' '{print $NF}'); ABRV_NAME=$(echo "$FILE" | awk -F '.' '{print $1 "_" $2 "_COMBINEDREPS"}'); fi
CELL_LINE=$(echo "$ABRV_NAME" | awk -F '_' '{print toupper($1) }')
CONDITION=$(echo "$ABRV_NAME" | awk -F '_' '{print toupper($2) }')
LOW_PVALUE="$ABRV_NAME""_p_1E-03"
DAY=$(echo $CONDITION | grep -Eo 'DAY[0-99]|DAY_[0-99]|D[0-99][^_]*');


#ASSIGN CONTROLS (INPUT)
if [ -z $INPUT ]; then echo "INPUT WAS NOT SET IN SUBMIT SCRIPT"
    if [[ "$CELL_LINE" = "MCF7" ]]; then
    INPUT="/users/j/a/jargordo/scratch/inputs/MCF7_genomic_input_TM_R1_HG19_33592775.bowtie"
    elif [[ "$CELL_LINE" = "MCF10a" ]]; then
        INPUT="/users/j/a/jargordo/scratch/inputs/MCF10A_genomic_input_TM_29M_reads.bowtie"
    elif [[ "$CELL_LINE" = "MDA" ]]; then
        INPUT="/users/j/a/jargordo/scratch/inputs/MDA231_genomic_input_TM_32M_reads.bowtie"
    elif [[ "$CELL_LINE" = "MDA231" ]]; then
        INPUT="/users/j/a/jargordo/scratch/inputs/MDA231_genomic_input_TM_32M_reads.bowtie"
    elif [[ "$CELL_LINE" = "BMSC" ]]; then
        if [[ "$DAY" = "DAY00" || "$DAY" = "D0" || "$DAY" = "D00" ]]; then
            INPUT="/users/j/a/jargordo/scratch/inputs/BMSC_DAY00_COMBINEDINPUT_48352613.bowtie"
        elif [[ "$DAY" = "DAY07" || "$DAY" = "D07" ]]; then
            INPUT="/users/j/a/jargordo/scratch/inputs/BMSC_DAY07_COMBINEDINPUT_27311093.bowtie"
        elif [[ "$DAY" = "DAY14" || "$DAY" = "D14" ]]; then
            INPUT="/users/j/a/jargordo/scratch/inputs/BMSC_DAY14_COMBINEDINPUT_54142153.bowtie"
        elif [[ "$DAY" = "DAY21" || "$DAY" = "D21" ]]; then
            INPUT="/users/j/a/jargordo/scratch/inputs/BMSC_DAY21_COMBINEDINPUT_47838179.bowtie"
        else echo "Can't set INPUT"
        fi
    else echo "Can't set INPUT"
    fi
else echo "INPUT WAS SET BY SUBMIT SCRIPT"
fi
echo "INPUT IS SET TO $INPUT"

#NOTE fix this for mouse submits

#ASSIGN OTHER VARIABLES
if [ -z $IDR ]; then echo "IDR WAS NOT SET IN SUBMIT SCRIPT. WILL GENERATE BOTH (LOW PVALUE AND NORMAL PVALUE) FILES"; IDR=2; else echo ""; fi
if [ -z $PVALUE ]; then echo "PVALUE WAS NOT SET IN SUBMIT SCRIPT. WILL USE DEFAULT (1e-05)"; PVALUE="1e-05"; else echo ""; fi
if [ $NOSPMR == 0 ]; then SPMR="-B --SPMR" else echo "BEDGRAPH PILEUP WAS DISABLED IN THE SUBMIT SCRIPT"; SPMR=""; fi
MFOLD="$MFOLDL $MFOLDU"
if [ -z $MFOLD ]; then MFOLD="5 50"; else echo "MFOLD WAS SET TO $MFOLD IN THE SUBMIT SCRIPT"; fi 

echo "------------------------MACS2_CALL_PEAKS_FOR_$ABRV_NAME-----------------------------------"
if [ $IDR == 0 ]; then
    macs2 callpeak -t $BOWTIE $BOWTIE2 -c $INPUT -m $MFOLD -f SAM -n $ABRV_NAME -g $GEN -p $PVALUE $SPMR
elif [ $IDR == 1 ]; then
    macs2 callpeak -t $BOWTIE -c $INPUT -m $MFOLD -f SAM -n $LOW_PVALUE -g $GEN -g $GEN -p 1e-03 --to-large
else
    macs2 callpeak -t $BOWTIE $BOWTIE2 -c $INPUT -m $MFOLD -f SAM -n $ABRV_NAME -g $GEN -p $PVALUE -B $SPMR
    macs2 callpeak -t $BOWTIE -c $INPUT -f SAM -m $MFOLD -n $LOW_PVALUE -g $GEN -g $GEN -p 1e-03 --to-large
fi


#6B)-----NORMALIZE_USING_FOLD_ENRICHMENT_AND_LOG_LR--------------------

#OPTION:-pseudocount (-p) to adjust for zero values is set to 0.00001

echo "------------------------MACS2_bdgcmp-----------------------------------"
MACS2_BEDGRAPH="$ABRV_NAME""_treat_pileup.bdg"
MACS2_CONTROL="$ABRV_NAME""_control_lambda.bdg"
rm $MACS2_CONTROL
#MACS2_OUTFILE2="$ABRV_NAME""_fold_enrichment.bdg"
#MACS2_OUTFILE3="$ABRV_NAME""_logLR.bdg"
#macs2 bdgcmp -t $MACS2_BEDGRAPH -c $MACS2_CONTROL -o $MACS2_OUTFILE1 -m FE -p 0.00001
#macs2 bdgcmp -t $MACS2_BEDGRAPH -c $MACS2_CONTROL -o $MACS2_OUTFILE2 -m logLR -p 0.00001

#7)-----CONVERT_BEDGRAPH_TO_BIGWIG-------------------------

#OPTION: generate random colour for wig tracks
#GRAPHCOLOUR=$((RANDOM %256 + 1)),$((RANDOM %256 +1)),$((RANDOM %256 +1))

#COLOR BREWER 12 COLOR QUALITATIVE PALETTE
C1="166,206,227";C2="31,120,180";C3="178,223,138";C4="51,160,44";C5="251,154,153";C6="227,26,28";C7="253,191,111";
C8="255,127,0";C9="202,178,214";C12="106,61,154";C11="255,255,153";C12="177,89,40"

if [ "$CONDITION" = "H3K4ME3" ]; then COLOUR=$C1;
    elif [ "$CONDITION" = "H3K4AC" ]; then COLOUR=$C2;
    elif [ "$CONDITION" = "H3K9ME3" ]; then COLOUR=$C3;
    elif [ "$CONDITION" = "H3K9AC" ]; then COLOUR=$C4;
    elif [ "$CONDITION" = "H3K27ME3" ]; then COLOUR=$C5;
    elif [ "$CONDITION" = "H3K27AC" ]; then COLOUR=$C6;
    elif [ "$CONDITION" = "H3K36ME3" ]; then COLOUR=$C7;
    elif [ "$CONDITION" = "H4K12AC" ]; then COLOUR=$C8;
    elif [ "$CONDITION" = "H4K20ME3" ]; then COLOUR=$C9;
    elif [ "$CONDITION" = "DAY00" ]; then COLOUR=$C1;
    elif [ "$CONDITION" = "DAY07" ]; then COLOUR=$C2;
    elif [ "$CONDITION" = "DAY14" ]; then COLOUR=$C3;
    elif [ "$CONDITION" = "DAY21" ]; then COLOUR=$C4;
else echo "$CONDITON not in the normal list: RANDOM COLOUR!"; COLOUR=$((RANDOM %256 + 1)),$((RANDOM %256 +1)),$((RANDOM %256 +1))
fi



if [ $NOSPMR == 0 ]; then
    for bedgraph in $MACS2_BEDGRAPH $MACS2_OUTFILE1; do
        FINAL_NAME=$(echo "$bedgraph" | awk -F '.' '{print $1}')
        CLIPPED="$FINAL_NAME"".clipped"
        FINALWIG="$FINAL_NAME"".bw"

        echo "GENERATING BIGWIG FOR $bedgraph"
        bedClip -verbose=2 $bedgraph $CHROMSIZE $CLIPPED
        bedGraphToBigWig $CLIPPED $CHROMSIZE $FINALWIG
        rm $CLIPPED

        gzip $bedgraph


#7A)-----GENERATE_CONTROL_FILE_FOR_UCSC_GENOME_BROWSER-------------------------

echo "track type=bigWig" "name=$FINAL_NAME" "description="'"'"MACS2 big wig from $FINAL_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/TM/$FINALWIG" "autoScale=off viewLimits=0.001:2 smoothingWindow=2 windowingFunction=maximum" "maxHeightPixels=128:48:8" "color=$COLOUR" "visibility=dense" >> UCSC_CONTROL_FILE_BIGWIGS_$DATE.txt
done
fi

#8)-----CONVERT_BEDS_TO_BIGBED-------------------------

PEAKSBED="$ABRV_NAME""_peaks.bed"
SUMMITBED="$ABRV_NAME""_summits.bed"
EXTENDSUMMIT="$ABRV_NAME""_extenedsummits.bed"
LPVPEAKSBED="$LOW_PVALUE""_peaks.bed"
#extend the summits
awk '{print $1 "\t" ($2 - 100) "\t" ($3 + 100) "\t" $4 "\t" $5}' $SUMMITBED > $EXTENDSUMMIT

for bedfile in $PEAKSBED $SUMMITBED $EXTENDSUMMIT $LPVPEAKSBED; do
CLIPPED="$PBS_JOBID"".clipped"
SORTED="$PBS_JOBID"".sorted"
FINAL_NAME=$(echo "$bedfile" | awk -F '.' '{print $1}')
FINALBED="$FINAL_NAME"".bb"
sort -k1,1 -k2,2n $bedfile > $SORTED
bedClip -verbose=2 $SORTED $CHROMSIZE $CLIPPED
bedToBigBed -type=bed3+2 $CLIPPED $CHROMSIZE $FINALBED
mv $CLIPPED $bedfile
rm $SORTED

#8A)-----GENERATE_CONTROL_FILE_FOR_UCSC_GENOME_BROWSER-------------------------
echo "track type=bigBed" "name=$FINAL_NAME" "description="'"'"peaks/summits from $FINAL_NAME"'"' "bigDataUrl=http://www.med.uvm.edu/steinlab/hg19/TM/$FINALBED" >> UCSC_CONTROL_FILE_BIGBED_$DATE.txt
done

echo "------------------------RUN_INFO----------------------------------"
echo "This is job ran on " `hostname` "on" `date +%d_%m_%Y_%H:%M:%S`



#Command to connect VACC to COMIS (you run this from vacc)
#smbclient \\\\med23.med.uvm.edu\\Steinlab-Research$ --user=med\\yourID