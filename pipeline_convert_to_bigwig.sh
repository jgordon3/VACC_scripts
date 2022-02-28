#!/bin/sh
#  pipeline_convert_to_bigwig.sh
#$ -cwd
#$ -j Y
# Request [num] amount of [type] nodes
#$ -pe threads 1


# var
TREAT=$1
GENOME="hg38" #conditional later
CHROMSIZE=$GENOME".chrom.sizes"
if [ ! -f $CHROMSIZE ]; then wget -np -nd -r http://genome.ucsc.edu/goldenpath/help/$CHROMSIZE; fi

# Naming and folders
ABRVNAME=${INPUT%.bdg}
CLIPPED=$ABRVNAME".clipped.bedgraph"
FINALNAME=$ABRVNAME".bw" #this will need to changed for multiple input file types
PUBLICFOLDER=/slipstream/galaxy/production/galaxy-dist/static/UCSCtracks/hg38/ChIP_seq




#samtools sort accepted_hits.bam sorted_hits
#genomeCoverageBed -bg -ibam sorted_hits.bam -g $CHROMSIZE > sorted_hits.bedgraph

bedClip -verbose=2 $INPUT $CHROMSIZE $CLIPPED
bedGraphToBigWig clipped.bedgraph $CHROMSIZE $FINALNAME


###### CREATE NORMALIZED BEDGRAPHS FROM BOWTIE SAM FILES##############
# these can be generated in MACS2 with --SPMR flag ... but I don't trust it
# count mapped reads
TREATMAPPEDREADS=$(samtools view -F 0x904 -c $TREAT);
CONTROLMAPPEDREADS=$(samtools view -F 0x904 -c $CONTROL)
printf "%s\n" "$DIVIDER" "MAPPED READS IN TREATMENT: $TREATMAPPEDREADS" "MAPPED READS IN CONTROL; $CONTROLMAPPEDREADS" "" "$DIVIDER"

# Scaling factor default is counts/million mapped reads
SCALE=1000000
TREATSCALEFACTOR=`echo "$SCALE/$TREATMAPPEDREADS" | bc -l`
CONTROLSCALEFACTOR=`echo "$SCALE/$CONTROLMAPPEDREADS" | bc -l`

# Normalize by scale factor
bedtools genomecov -ibam $TREAT -bg -scale $TREATSCALEFACTOR -g $CHROMSIZE >> $TREATBDG
bedtools genomecov -ibam $CONTROL -bg -scale $CONTROLSCALEFACTOR -g $CHROMSIZE >> $CTRLBDG

# Convert bedgraphs to bigwig
wigToBigWig -clip $TREATBDG $CHROMSIZE $TREATWIG
wigToBigWig -clip $CTRLBDG $CHROMSIZE $CTRLWIG
wigToBigWig -clip $TREATFE $CHROMSIZE $TREATFEWIG
wigToBigWig -clip $TREATlogLR $CHROMSIZE $TREATlogLRWIG



#Cleanup
rm clipped.bedgraph

cp $FINALNAME $PUBLICFOLDER



https://galaxy.med.uvm.edu/static/UCSCtracks/hg38/ChIP_seq/Runx2.imsc3.diff.MACS2_in_Galaxy_treat_pileup.bw
