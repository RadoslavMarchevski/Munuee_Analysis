#!/bin/bash

################################################
#  This is a simple bash script which submites #
#  jobs on Batch                               #
################################################

type="data";
ver=1
queue="short"

dir="/home/user/compact/reader/LSF/$ver/$type/";
echo $dir
index=0;
i=1;
ls "$dir" | while read -r file; 
  do 
  bfile="$file"
  echo $bfile;
  cd $dir
  #bsub -q $queue -o /dev/null -J "$type.$ver.$index" $file
  bsub  -W 300 -app Reserve300M -n 1 -q $queue -o /gpfs/fs3/na62/user/user/LSF/$ver/$type/$type.$ver.$index.out -J "./$type.$ver.$index" $dir/$file
  echo "bsub  -W 180 -app Reserve300M -n 1 -q $queue -o /gpfs/fs3/na62/user/user/LSF/$ver/$type/$type.$ver.$index.out -J "$type.$ver.$index" $dir/$file"
  ((index+=i))
done

