#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=14gb,pvmem=15gb
#PBS -l walltime=30:00:00
#PBS -N BamtoBigwig
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea
if [ -z $PARENTDIR ]; then echo "NO DIRECTORY SET"; exit 1; else echo ""; fi
if [ -z $CHROMSIZE ]; then echo "NO CHROMESIZE SET"; exit 1; else echo ""; fi

cd $PARENTDIR
FINAL="$FOLDER""_accepted_hits.bw"
cd $FOLDER

samtools sort accepted_hits.bam sorted_hits
genomeCoverageBed -bg -ibam sorted_hits.bam -g $CHROMSIZE > sorted_hits.bedgraph
bedClip -verbose=2 sorted_hits.bedgraph $CHROMSIZE clipped.bedgraph
bedGraphToBigWig clipped.bedgraph $CHROMSIZE $FINAL
rm clipped.bedgraph sorted_hits.bedgraph

cd $PARENTDIR
