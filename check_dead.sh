#!/bin/bash

tag=$1
mkdir test/$tag
for i in `ls | grep $tag` ;
do
    #echo $PWD
    
    cat $i/compact.out | tail -4 >> test/$tag/check_$tag.out
    error=`cat $i/compact.out | grep 'SCMP' | wc -l`
    #echo $error
    #if [ $error -eq 1 ] ; then
    #	#echo "cat $i/compact.out | grep 'bad_alloc'"
    #	cat $i/compact.out | grep 'bad_alloc' > test/$tag/check_good_${tag}_${i}.out; 
    #fi
    
    #if [ $error -eq 0 ] ; then
    #	#echo "cat $i/compact.out | grep 'bad_alloc'"
    #	cat $i/compact.list  > test/list_fail_${tag}_${i}.list; 
    #	#rm /afs/cern.ch/user/r/rmarchev/compact/reader/$i/run2007.root
    #fi
done
