#!/bin/sh
#run a single cell expeiment on cellranger (10x genomics)
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/mkfastq
#installed locally at: /slipstream/home/jonathan/bin/cellranger-3.0.2/cellranger
alias cellranger='./slipstream/home/jonathan/bin/cellranger-3.0.2/cellranger'
#SEE usage for required arguments:
#

usage ()
{
case $ERROR in
    1) MESSAGE="You did not supply a valid directory with the --run flag or it does not exist";;
    2) MESSAGE="You did not supply a valid csv with the --csv flag or it does not exist";;
    *) MESSAGE="Help file activated or Random Error: resubmit";;
esac
printf "%s\n" "$DIVIDER" "" "$MESSAGE" "" "This is script is now submitted via qsub with the command" "";
printf "%s\n" "qsub run_cellranger_on_SGE.sh --run /path/to/fastq" ""
printf "%s\n" "****  optional flags: ******** " "" "$DIVIDER" ""
exit 1; #finish this later

id=$1   #random name for the experiment (sets up folder in dir from where the command is called)
run=$2  #location of stock Illumina bcl folder
csv=$3  #simple csv samplesheet

cellranger mkfastq --id=$id --run=$run --csv=$csv --qc



