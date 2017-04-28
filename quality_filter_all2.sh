#!/bin/bash
filelist=($(ls *.fastq))
for i in "${filelist[@]}"
do 
   if [ ! -f filtered/$i ];
   then
       fastq_quality_filter -q 24 -p 90 -i $i -o filtered/$i
   fi
done
