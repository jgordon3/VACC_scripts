#PBS -l nodes=1
#PBS -l walltime=01:03:00
#PBS -joe

cd $PBS_O_WORKDIR
NOW=$(date +"%m/%d/%y %H-%M")
echo $NOW

# launch continuation job

echo "starting jobtwo job" >> count_to_num.out

# Again, redirect output to append to count_to_num.out, put
# job in background
cr_restart jobone.checkpoint >> count_to_num.out 2>&1 &

BGPID=$!

# this number should be adjusted for longer jobs
echo "Sleeping 140 seconds to allow to run..."
sleep 1800

NOW=$(date +"%m/%d/%y %H-%M")
echo $NOW

cr_checkpoint -p $BGPID -f jobone.checkpoint --term