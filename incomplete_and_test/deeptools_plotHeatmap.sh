#!/bin/sh

#  deeptools_plotHeatmap.sh
#  
BIGWIG=$1
BIGWIG=$(readlink -f $BIGWIG)
if [ ! -f "$BIGWIG" ]; then echo "Need valid path to bigwig file"; fi
BASE=`basename $BIGWIG`

ABRV_NAME=${BASE%.fastq}

FILT_FASTQ="${ABRV_NAME}_trim_filt.fastq";
bed=$2


ABRVNAME=$(echo $bigwig | awk -F '_' '{print $1"_"$2}')
echo $ABRVNAME
compmat="$ABRVNAME.computematrix"
plotout="$ABRVNAME.png"

BEDNAME=$(echo $bed | awk -F '_' '{print $1"_"$2}')
reglabel="$BEDNAME_peaks"

computeMatrix reference-point -S $bigwig -R $bed -a 1000 -b 1000 -out $compmat --referencePoint center


plotHeatmap -m $compmat -out $plotout --heatmapHeight 15 --refPointLabel peak.center --regionsLabel "$reglabel" --plotTitle "$ABRVNAME"

