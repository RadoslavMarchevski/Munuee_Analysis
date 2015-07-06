
export inputtag=$1
export tag=$2
export outversion=$3

#export outlogdir=/afs/cern.ch/user/r/rmarchev/compact/reader/test/$inputtag
export outputdir=/afs/cern.ch/user/r/rmarchev/compact/reader/test/$inputtag
#export outputdirall=$outputdir/${tag}_$outversion/all

echo "source /home/rmarchev/compact/reader/merge_exrun_1.sh" | bsub -W 300 -app Reserve3G -q short  -J ${inputtag}_${tag}_${outversion}  -o $outlogdir/${tag}_$outversion/merged_${tag}_${outversion}.out -e $outlogdir/${tag}_$outversion/merged_${tag}_${outversion}.err


