#!/usr/bin/env python

#  HTSeq.py
#  
#
#  Created by Jonathan Gordon on 4/23/14.
#

import HTSeq
features =  HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for line in open( "DAY14_COMBREPS_42837018_peaks_IDR_filtered.htseqin" ):
    fields = line.split( "\t" )
    iv = HTSeq.GenomicInterval( fields[1], int(fields[2]), int(fields[3]) )
    features[ iv ] += fields[0]

import collections
counts = collections.Counter( )

almnt_file = HTSeq.SAM_Reader( "BMSC_RUNX2D14_REP2_18977368.bowtie" )
for almnt in almnt_file:
    if not almnt.aligned:
        counts[ "_unmapped" ] += 1
        continue
    gene_ids = set()
    for iv, val in features[ almnt.iv ].steps():
        gene_ids |= val
    if len(gene_ids) == 1:
        gene_id = list(gene_ids)[0]
        counts[ gene_id ] += 1
    elif len(gene_ids) == 0:
        counts[ "_no_feature" ] += 1
    else:
        counts[ "_ambiguous" ] += 1

for gene_id in counts:
    print gene_id, counts[ gene_id ]
#f = open('HTSeq.out','w')
#f.write(gene_id)
#f.write("\t")
#f.write(counts[ gene_id ])
#f.close()

