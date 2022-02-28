#!/bin/bash
#$ -q slipstream_queue@cn-t630

#use: move this script and a samplesheet (script will search for file with the suffix .csv) to Illumina run directory (YYMMDD_machinename_NNNN). 
#submit script with path to run directory
#qsub run_demux_v3.sh PATH/TO/DIR  

DIR=$1
cd $DIR
INT_DIR="${DIR}/Data/Intensities"
DEMUX_REPORTS="${DIR}_demultiplex_reports"
FASTQ_DIR="./${DIR}_fastqs"
DEMUX_RESULT_DIR="${DIR}_DATA"


#find samplesheet
SAMP_SHEET=$(find $DIR -maxdepth 1 -name '*.csv')

#checks
if [ -f $SAMP_SHEET ]; then echo "Found samplesheet ==> looking for data"; else echo "Could not find SampleSheet.csv"; exit 1; fi
if [ -d $INT_DIR ]; then echo "Found data ==> starting configuration"; else echo "This folder does not appear to have Data to convert. Try relaunching script from a folder containing data"; exit 1; fi

#options (to flag later)
#barcode mismatch value
BC_MISMATCH=0

/slipstream/home/jonathan/bin/bcl2fastq2-v2.18.0.12/bin/bcl2fastq -p 12 --barcode-mismatches $BC_MISMATCH --sample-sheet $SAMP_SHEET --reports-dir $DEMUX_REPORTS --output-dir $FASTQ_DIR


#move stuff for easier extraction
mkdir "../${DEMUX_RESULT_DIR}"
mv $DEMUX_REPORTS "../${DEMUX_RESULT_DIR}/"
mv $FASTQ_DIR "../${DEMUX_RESULT_DIR}/"

#make pdf reports
#/usr/local/bin/wkhtmltopdf --orientation landscape index.html test.pdf

#compress for archiving
ARCH_NAME="${DIR}_PROCESSED.tar.gz"
tar -cvzf $ARCH_NAME $DIR && rm -R $DIR


