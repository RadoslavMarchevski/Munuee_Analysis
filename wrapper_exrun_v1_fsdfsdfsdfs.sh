outversion=$1
nofevents=$2
queue=$3
tag=$4

if [ ! -z "$5" ] ; then lcount=$5;	else	lcount=0;	fi
if [ ! -z "$6" ] ; then hcount=$6;	else	hcount=0;	fi

noffilesperrun=10

#inputtag=ee_split_2004
#inputdir=/gpfs/fs3/na62/na48/data/2004/$inputtag
inputtag=ee_split_2003
#inputtag=ee_pass5
inputdir=/gpfs/fs3/na62/na48/data/2003/$inputtag #ask gia !!!
#inputdir=/gpfs/fs3/na62/na48/data/2003/ee_split_2003/
#inputtag=k2pigd_kloe_pi0d_radcor
#inputdir=/gpfs/fs3/na62/na48/mc/cmc007_v2/$inputtag

export inputlistsdir=/home/rmarchev/compact/reader/tmp/lists/$inputtag
export outlogdir=/home/rmarchev/compact/reader/tmp/logs/$inputtag
export outputdir=/home/rmarchev/compact/reader/tmp/output/$inputtag
export outputdirall=/home/rmarchev/compact/reader/tmp/output/$inputtag/${tag}_$outversion/all

if [ ! -d $inputlistsdir ]; 	then mkdir $inputlistsdir; 	fi 
if [ ! -d $outlogdir ]; 	then mkdir $outlogdir; 		fi 
if [ ! -d $outputdir ]; 	then mkdir $outputdir; 		fi 

if [ ! -d $inputlistsdir/${tag}_$outversion ]; 	then mkdir $inputlistsdir/${tag}_$outversion; 	fi 
if [ ! -d $outlogdir/${tag}_$outversion ]; 	then mkdir $outlogdir/${tag}_$outversion; 	fi 
if [ ! -d $outputdir/${tag}_$outversion ]; 	then mkdir $outputdir/${tag}_$outversion; 	fi 
if [ ! -d $outputdir/${tag}_$outversion/all ]; 	then mkdir $outputdir/${tag}_$outversion/all; 	fi 


noffiles=`ls $inputdir | wc -l`
let "nofruns=$noffiles/$noffilesperrun + 1"

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

		ls -1 $inputdir/*  | head -$fihead | tail -$fitail > $inputlistsdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.list
		
		if [ ! -d $outputdir/${tag}_$outversion/${i} ];       then mkdir $outputdir/${tag}_$outversion/${i};        fi
	
		#bsub -W 300 -app Reserve1G -q short -o /home/rmarchev/compact/reader/rado.txt /home/rmarchev/compact/reader/simple.sh 
	

		echo "source $HOME/compact/reader/simple.sh $i $outversion $nofevents $tag" | bsub -W 300 -app Reserve1G -q $queue  -J ${inputtag}_${i}  -o $outlogdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.out -e $outlogdir/${tag}_$outversion/${tag}_${outversion}_run_${i}.err	
	done




