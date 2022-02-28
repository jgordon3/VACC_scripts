#!/bin/bash
#PBS -l nodes=1
#PBS -l walltime=00:03:00
#PBS -j oe

echo $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

# launch job
echo "starting jobone to output file" >> count_to_num.out

# Note that I have redirected both STDERR and STDOUT to count_to_num.out.  This is to make
# sure that these output streams go somewhere that cr_checkpoint can checkpoint.  If you don't
# redirect, you will discover that you can't restart the job on a different node.
# Also, note that it's running in the background, allowing the script to continue while it runs.

cr_run $HOME/bin/scripts/submit_CUFFDIFF_ISHIKAWA_all.sh >> count_to_num.out 2>&1 &

BGPID=$!
echo $BGPID
# Sleep for most of, but not all of, the walltime.  You want to allow your job to have time to
# checkpoint before the walltime is up.  For a large job with a large amount of memory, you
# might want to allow a long time to checkpoint.
echo "Sleeping 140 seconds to allow to run..."
sleep 140

# This will checkpoint, then terminate the job
cr_checkpoint -p $BGPID -f jobone.checkpoint --term