#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=8gb,pvmem=9gb
#PBS -l walltime=24:00:00
#PBS -N gunzip
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

FILEPATH=/users/j/a/jargordo/scratch/sequences/
FILES=/users/j/a/jargordo//scratch/sequences/*.gz
cd $FILEPATH

for f in $FILES; do
gzip $f

done



