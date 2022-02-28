#!/bin/bash
#PBS -l nodes=1:ppn=1,pmem=1gb,pvmem=1gb
#PBS -l walltime=00:10:00
#PBS -N TEST
#PBS -j oe
#PBS -M Jonathan.A.Gordon@uvm.edu
#PBS -m bea

echo "this worked"
echo "$VAR"
