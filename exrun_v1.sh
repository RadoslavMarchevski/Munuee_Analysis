#LSF=.

runid=$1
outversion=$2
nofevents=$3
tag=$4


#source /home/rmarchev/compact/reader/var.sh

if [ $nofevents -gt 0 ] ; then
    /afs/cern.ch/user/r/rmarchev/compact/reader/./compact -l $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${runid}.list  -rootfilename /${tag}_${outversion}_run_${runid} -nevt $nofevents #-cheat
else
    /afs/cern.ch/user/r/rmarchev/compact/reader/./compact -l $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${runid}.list  -rootfilename ${tag}_${outversion}_run_${runid} #-cheat
fi

cp -f ${tag}_${outversion}_run_${runid}.root  $outputdir/${tag}_$outversion/${runid}/${tag}_${outversion}_run_${runid}.root
rm  ${tag}_${outversion}_run_${runid}.root





