#!/bin/sh

#finds replicate pairs


for f in MCF7*; do
    CELL=$(echo $f | awk -F '_' '{print $1}')
    HISTONE=$(echo $f | awk -F "_" '{print $2}')
    R1=$CELL"_"$HISTONE"_CONTROL_R1"*
        if [ -e $R1 ]; then
            R2=$CELL"_"$HISTONE"_CONTROL_R2"*
            CHECK="$R1""_""$R2"
                if [[ "$CHECK" = "$CHECK2" ]];
                    then echo "Already done"
                    else echo "send $R1 $R2"; fi
            else echo "no matching replicate"; fi
CHECK2=$CHECK
done



