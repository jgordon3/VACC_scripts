#!/bin/sh

#  run_replicate_analysis_on_SGE.sh
#  
#
#  Created by Jonathan Gordon on 11/13/15.
#

#bamCorrelate

bamCorrelate bins --bamfiles MCF10A_H3K27AC_R1_06202013_sorted.bam MCF10A_H3K27AC_R2_03112014_sorted.bam --labels MCF10A_R1 MCF10A_R2 --binSize 1000 --corMethod spearman -f 200 --colorMap Reds --zMin 0.5 --zMax 1 -o correlation_spearman.pdf