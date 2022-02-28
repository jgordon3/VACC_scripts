#!/bin/sh

#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

# Pass variables from MACS2 script

# OPTIONS
while getopts "T:I:cghtusa" opt; do
    case $opt in
        T) TREAT=$OPTARG;;          #fastq is required
        I) CONTROL=$OPTARG;;
        c) COMPRESS=$OPTARG;;
        g) GENOME=$OPTARG;;
        p) PVALUE=$OPTARG;;
        u) NAME=$OPTARG;;
        s) SAVEBED=$OPTARG;;
        a) ADD_PEAKCALL=$OPTARG;;
        h) HELP=1;;
        :) HELP=1;;
        *) HELP=1;;
        esac
done
if [[ -e $TREAT ]]; then TREAT=$TREAT; else echo "Need valid path to BOWTIE BAM/SAM file"; HELP=1; fi
if [[ -e $CONTROL ]]; then CONTROL=$CONTROL; else echo "No Input set"; HELP=1; fi
if [[ ! -z "$COMPRESS" ]]; then COMPRESS=$COMPRESS; else COMPRESS="y"; fi
if [[ ! -z "$GENOME" ]]; then GENOME=$GENOME; else GENOME="hg38"; fi
if [[ ! -z "$PVALUE" ]]; then PVALUE=$PVALUE; else PVALUE="1e-05"; fi
if [[ ! -z "$NAME" ]]; then NAME=$NAME; else NAME="J. Gordon"; fi
if [[ ! -z "$SAVEBED" ]]; then SAVEBED=$SAVEBED; else SAVEBED="y"; fi #save bedgraphs option
if [[ ! -z "$ADD_PEAKCALL" ]]; then ADD_PEAKCALL=$ADD_PEAKCALL; else ADD_PEAKCALL="n"; fi
if [[ -z $HELP ]]; then echo "";
    else
        printf "%s\n" "" "This is script is now submitted via qsub with the command" "" "qsub run_MACS2_on_SGE.sh -T /path/to/Treatment/BAM -I /path/to/input/control/BAM" "" "****  optional flags: ******** " "" "-c: [y/n] Compress outputs as gzip/tar archives. Default is y" "-g: Genome for alignment. Default is hg38. Some other options are hg19, mm9, mm10 ... etc" "-p: Set p value. Default is 1e-05" "-u: User info for README file" "-s: [y/n] Save bedgraph outputs. Bigwigs are generated, so deleting redundant bedgraphs saves some space. Default is Y" "-a: [y/n] Launch additional broadpeak calls using MACS2 broadpeak and SICER. Default in n" ""
    exit 1;
fi
shift $(($OPTIND -1))

#ASSIGN OTHER VARIABLE ---- put some conditionals in here
TREATABRVNAME=${TREAT%.bam}
CONTROLABRVNAME=${CONTROL%.bam}
SICER_BIN=/slipstream/home/jonathan/bin/SICER_V1.1


mkdir $TREATABRVNAME
cd $TREATABRVNAME
cp -r $SICER_BIN ./

sh DIR/SICER.sh ["InputDir"] ["bed file"] ["control file"] ["OutputDir"] ["Species"] ["redundancy threshold"] ["window size (bp)"] ["fragment size"] ["effective genome fraction"] ["gap size (bp)"] ["FDR"]


