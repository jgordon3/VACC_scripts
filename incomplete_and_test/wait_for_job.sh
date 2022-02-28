#!/bin/bash
# wait_for_job.sh
# waits for a qsub job to complete
STARTTIME=$(date +%s)

# Matches a list of JOBIDs to currently running /queued jobs in SGE
qstat_query () {
    SUBMITS="$1" # Note: it is not easy to pass an array to a function in bash so the JOBIDS are passed as a string variable
    unset SUBMIT_ARRAY # clear previous
        for i in $SUBMITS; do SUBMIT_ARRAY+=($i); done
        echo "Submitted: ${#SUBMIT_ARRAY[@]} jobs: ${SUBMIT_ARRAY[@]}"

    JOB_MATCH=1
    while [[ $JOB_MATCH -gt 0 ]]; do
        sleep 15    #change interval
        QSTAT_QUERY=`qstat -u "*"`
        JOB_QUERY=$(echo $QSTAT_QUERY |grep -Eo '[0-9]{5}') # not optimal .. returns all 5 digit numbers in qstat

        unset JOBIDS_QUERY # clear previous query and list running jobs
            for i in $JOB_QUERY; do JOBIDS_QUERY+=($i); done

        JOB_MATCH=0  #match jobs
            for i in "${SUBMIT_ARRAY[@]}"; do
                for j in "${JOBIDS_QUERY[@]}"; do
                    if [[ $i == $j ]]; then JOB_MATCH=$((JOB_MATCH+1));fi
                done
            done
        CURRENTTIME=$(date +%s); ELAPSED=$(($CURRENTTIME - $STARTTIME)); echo "checking for job(s): ${SUBMIT_ARRAY[@]}"; echo "$JOB_MATCH are still running. $ELAPSED seconds elapsed"
    done
echo "all jobs complete"
}


# LAUNCH QSUB for JOB and catch new JOBID
# Could make this a function if needed

SCRIPT=wait.sh
FILELIST=""
while IFS='' read -r line || [[ -n "$line" ]]; do FILELIST="$FILELIST $line"; done < test.txt
#cut columns?

# submit all jobs and return a list of jobs
for f in $FILELIST; do
    JOB1=$(qsub $SCRIPT $f) #add file here
    #returns whole line
    JOBID=$(awk -v RS=[0-9]+ '{print RT+0;exit}' <<< "$JOB1") #returns JOBID
    JOBIDS="$JOBIDS $JOBID"
done


qstat_query "$JOBIDS"

rm wait.sh.*


