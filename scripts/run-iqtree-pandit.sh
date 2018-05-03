#! /bin/bash

for ((i=1;i<10000;i++)); do
        if [ ! -d $i ]; then
                continue
        fi
        cd $i
        aln=data.$i
        model=`cat model.$i`
        echo $aln $model
        submit2sge -q cluster -r old -N sbs$i "iqtree -s $aln -m $model -b 100"
        cd ..
done
