# merges fastqs from individual folders
# integrate into demux script

parentpath=`pwd`
parentfolder=`basename $parentpath`
finalfolder="$parentfolder""_merged_fastqs"
finalpath=$parentpath/$finalfolder
echo $finalpath
mkdir $finalpath

for f in *; do    #folders
    if [ -d "$f" ]; then
        cd $parentpath/$f
        fn=$(ls | head -n1 | awk -F "." '{print $1}'); # list of files
        final=$(echo $fn | awk -F "_00*" '{print $1"_cat000.fastq.gz"}');
            for x in *.fastq.gz; do
            echo "$x copied to $final" #>> $finalpath/MERGE_LOG.txt
            echo $final
            #cat $x >> $final
            done
#        cp $final $finalpath
        cd $parentpath
    fi
done


# lane specific

for f in *; do #folder
    L001=`ls $f/*_L001_*`;
    L002=`ls $f/*_L002_*`;
    for i in $L001; do
        base=${i%_0**.fastq.gz}
        base="$base""_cat000.fastq.gz"
        cat $i >> $base
        echo "$i catted to $base" >> cat_log.txt
    done;
    for i in $L002; do
        base=${i%_0**.fastq.gz}
        base="$base""_cat000.fastq.gz"
        cat $i >> $base
        echo "$i catted to $base" >> cat_log.txt
    done;
done
