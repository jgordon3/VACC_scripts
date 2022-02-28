 #! bin/bash

 alias bcl2fastq="~/bin/bcl2fastq2-v2.18.0.12/bin/bcl2fastq"  

 #FLOWCELL_DIR=$1
FLOWCELL_DIR=/slipstream/home/jonathan/160919_SNL128_0144_AH3GMLBCXY
OUTPUT_DIR=~/TEST
SAMPLE_SHEET_PATH=~/Samplesheet_bcl2fastq2_format_09272016.csv

 bcl2fastq --use-bases-mask=Y150,I8,Y150 \
  --create-fastq-for-index-reads \
  --minimum-trimmed-read-length=8 \
  --mask-short-adapter-reads=8 \
  --ignore-missing-positions \
  --ignore-missing-controls \
  --ignore-missing-filter \
  --ignore-missing-bcls \
  -r 6 -w 6 \
  -R ${FLOWCELL_DIR} \
  --output-dir=${OUTPUT_DIR} \
  --interop-dir=${INTEROP_DIR} \
  --sample-sheet=${SAMPLE_SHEET_PATH}

   bcl2fastq -R ${FLOWCELL_DIR} --sample-sheet=${SAMPLE_SHEET_PATH} --output-dir=${OUTPUT_DIR}