
LSF=/jobdir/${LSB_JOBID}

runid=$1
outversion=$2
nofevents=$3
tag=$4
inputlistdir=$HOME/compact/reader/tmp/lists/
#tag=ee_pass5
outputdir=$HOME/compact/reader/tmp/output/

source /home/rmarchev/compact/reader/var.sh

if [ $nofevents -gt 0 ] ; then
	/home/rmarchev/compact/reader/./compact -l $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${runid}.list  -rootfilename ${LSF}/${tag}_${outversion}_run_${runid} -nevt $nofevents
else
	/home/rmarchev/compact/reader/./compact -l $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${runid}.list  -rootfilename ${LSF}/${tag}_${outversion}_run_${runid}
fi

cp -f ${LSF}/${tag}_${outversion}_run_${runid}.root  $outputdir/${tag}_$outversion/${runid}/${tag}_${outversion}_run_${runid}.root
rm  ${LSF}/${tag}_${outversion}_run_${runid}.root





