#!/bin/sh

#  Pipeline_metrics.sh
#  
#
#  Created by Jonathan Gordon on 2/22/17.
#

#text log file ouput from STAR
STAR_logfile=$1

#initial clean with sed

cat $STAR_logfile | sed -e 's/^[ \t]*//' | sed '/^\s*$/d'

