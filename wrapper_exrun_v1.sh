outversion=$1
nofevents=$2
queue=$3
tag=$4


if [ ! -z "$6" ] ; then lcount=$6;	else	lcount=0;	fi
if [ ! -z "$7" ] ; then hcount=$7;	else	hcount=0;	fi

noffilesperrun=2

#inputtag=ee_split_2004
#inputdir=/gpfs/fs3/na62/na48/data/2004/$inputtag
#inputtag=ee_split_2003
inputtag=pass5
#inputdir=/gpfs/fs3/na62/na48/data/2003/$inputtag #ask gia !!!
#inputlist=/afs/cern.ch/user/r/rmarchev/compact/reader/SC_SS0123_pass5_tot.list #ask
inputlist=/afs/cern.ch/user/r/rmarchev/compact/reader/try.list #ask
#inputdir=/gpfs/fs3/na62/na48/data/2003/ee_split_2003/
#inputtag=k2pigd_kloe_pi0d_radcor
#inputdir=/gpfs/fs3/na62/na48/mc/cmc007_v2/$inputtag

export inputlistsdir=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/lists/$inputtag
export outlogdir=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/logs/$inputtag
export outputdir=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/output/$inputtag
export outputdirall=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/output/$inputtag/${tag}_$outversion/all

if [ ! -d $inputlistsdir ]; 	then mkdir $inputlistsdir; 	fi
if [ ! -d $outlogdir ]; 	then mkdir $outlogdir; 		fi
if [ ! -d $outputdir ]; 	then mkdir $outputdir; 		fi

#mkdir $inputlistsdir;
#mkdir $outlogdir;
#mkdir $outputdir;
if [ ! -d $inputlistsdir/${tag}_$outversion ]; 	then mkdir $inputlistsdir/${tag}_$outversion; 	fi
if [ ! -d $outlogdir/${tag}_$outversion ]; 	then mkdir $outlogdir/${tag}_$outversion; 	fi
if [ ! -d $outputdir/${tag}_$outversion ]; 	then mkdir $outputdir/${tag}_$outversion; 	fi
if [ ! -d $outputdir/${tag}_$outversion/all ]; 	then mkdir $outputdir/${tag}_$outversion/all; 	fi


noffiles=`cat $inputlist | wc -l`
let "nofruns=$noffiles/$noffilesperrun"

if [ $lcount -eq 0 ] ; then lcount=1; 		fi
if [ $hcount -eq 0 ] ; then hcount=$nofruns;	fi

for i in `seq $lcount $hcount`;
	do
		let "fihead=$noffilesperrun*$i"

		if [ $i -eq $nofruns ];	then
			let "fitail=$noffiles-$noffilesperrun*($i-1)"
		else
			let "fitail=$noffilesperrun"
		fi

		cat $inputlist  | head -$fihead | tail -$fitail > $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.list

		if [ ! -d $outputdir/${tag}_$outversion/${i} ];       then mkdir $outputdir/${tag}_$outversion/${i};        fi



		#run bapp command to see available profiles
		echo "source /afs/cern.ch/user/r/rmarchev/compact/reader/exrun_v1.sh $i $outversion $nofevents $tag" | bsub -W 300  -q $queue  -J ${inputtag}_${i}  -M 3000000 -o $outlogdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.out -e $outlogdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.err


		#link the outpur in advance. Yhose runs, which failed to produce output are easy to find by dead links
		ln -sf $outputdir/${tag}_$outversion/${i}/${tag}_${outversion}_run_${i}.root  $outputdirall/${tag}_${outversion}_run_${i}.root
	done
