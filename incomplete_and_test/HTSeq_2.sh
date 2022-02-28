#!/usr/bin/env python

#runs HTSeq

os.getenv('PEAKFILE')
os.getenv('SAMFILE')
import HTSeq
features =  HTSeq.GenomicArrayOfSets( "auto", stranded=False )
for line in open( "PEAKFILE" ):
    fields = line.split( "\t" )
    iv = HTSeq.GenomicInterval( fields[1], int(fields[2]), int(fields[3]) )
    features[ iv ] += fields[0]

import collections
counts = collections.Counter( )

almnt_file = HTSeq.SAM_Reader( "SAMFILE" )
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
#print gene_id, counts[ gene_id ]
    f = open('HTSeq.out','w')
    f.write(gene_id)
    f.write("\t")
    f.write(str(counts[ gene_id ]))
    f.close()

