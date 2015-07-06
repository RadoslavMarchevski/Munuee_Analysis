
LSF=/jobdir/${LSB_JOBID}

mergedfile=merged_${inputtag}_${tag}_${outversion}_onbatch_2.root

infiles=""

for file in `ls $outputdirall`;
	do
		#ls $outputdirall/$file
		if [ -e $outputdirall/$file ]; then
	     		infiles+="$outputdirall/$file "
		fi
	done

echo $infiles
echo $mergedfile

cd $LSF
hadd $mergedfile $infiles

ls -l
cp -f $mergedfile $outputdir/${tag}_$outversion/$mergedfile
rm $mergedfile

