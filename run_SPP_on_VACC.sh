#!/bin/bash

#PBS -N SPP-[R-2.11.1]
#PBS -m be
#PBS -k oe

#USE R v2.11.1
/users/j/a/jargordo/bin/R-2.11.1/bin/R

# load the library
library(spp);

# The following section shows how to initialize a cluster of 8 nodes for parallel processing
# see "snow" package manual for details.
library(snow)
cluster <- makeCluster(8);

chip.data <- read.bowtie.tags("chip.bowtie.file");
input.data <- read.bowtie.tags("input.bowtie.file",max.bowtie.tag.length=32);

# get binding info from cross-correlation profile
# srange gives the possible range for the size of the protected region;
# srange should be higher than tag length; making the upper boundary too high will increase calculation time
#
# bin - bin tags within the specified number of basepairs to speed up calculation;
# increasing bin size decreases the accuracy of the determined parameters
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5,cluster=cluster);

R CMD INSTALL mypkg -l /my/own/R-packages/