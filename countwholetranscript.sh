#!/bin/bash
filelist1=(11.bam.sorted.bam
12.bam.sorted.bam
13.bam.sorted.bam
14.bam.sorted.bam
15.bam.sorted.bam
16.bam.sorted.bam
17.bam.sorted.bam
18.bam.sorted.bam
19.bam.sorted.bam
20.bam.sorted.bam
21.bam.sorted.bam
22.bam.sorted.bam
23.bam.sorted.bam
24.bam.sorted.bam
2.bam.sorted.bam
3.bam.sorted.bam
4.bam.sorted.bam
5.bam.sorted.bam
6.bam.sorted.bam
7.bam.sorted.bam
8.bam.sorted.bam
9.bam.sorted.bam
)
filelist2=(P408.bam.sorted.bam
P411.bam.sorted.bam
P414.bam.sorted.bam
P422.bam.sorted.bam
P426.bam.sorted.bam
901.bam.sorted.bam
904.bam.sorted.bam
910.bam.sorted.bam
913.bam.sorted.bam
918.bam.sorted.bam
)

for i in "${filelist1[@]}"
do
  echo $i
  samtools index $i
  samtools idxstats $i > wholetranscript/$i.out
  
done
for i in "${filelist2[@]}"
do
  echo $i
  samtools index $i
  samtools idxstats $i > wholetranscript/$i.out
done