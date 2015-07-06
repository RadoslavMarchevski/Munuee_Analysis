

dryrun=$1

outversion=v1
tag=K2pi	
inputtag=k2pigd_kloe_pi0d_radcor
#inputtag=ee_split_2004
#inputtag=ee_split_2003
#inputtag=ee_pass5
outputdirall=/home/rmarchev/compact/reader/tmp/output/$inputtag/${tag}_$outversion/all

lnum=0
hnum=0
submit=0
submit2=0

for i in $outputdirall/*; 
	do 
		if [ ! -s $i ]; then 
			file=${i/$outputdirall\//}; 
			file=${file/${tag}_${outversion}_run_/};
			num=${file/.root/}
			echo $num

			if [ $lnum -eq 0 ]; then
				lnum=$num
				hnum=$num
				continue
			fi 

			if [ $num -ge $lnum ]; then 
				hnum=$num
			else 
				submit=1
				submit2=1
			fi 
		else 
			submit=1
		fi

		if [ $submit -eq 1 ] ; then	
			
			if [ $lnum -ne 0 ]; then 
				echo $lnum " -- " $hnum
				
				if [ $dryrun -eq 0 ] ; then 
					./wrapper_exrun_v1.sh $outversion -1 long $tag $lnum $hnum  		
				fi

				if [ $submit2 -eq 1 ]; then	
					echo $num " -- " $num
					
					if [ $dryrun -eq 0 ] ; then 
						./wrapper_exrun_v1.sh $outversion -1 long $tag $num $num  		
					fi
				fi 	
			fi
	
			submit2=0	
			submit=0	
			lnum=0
			hnum=0
		fi
	done


if [ $lnum -ne 0 ]; then 
	echo $lnum " -- " $hnum
	
	if [ $dryrun -eq 0 ] ; then 
		./wrapper_exrun_v1.sh $outversion -1 long $tag $lnum $hnum  		
	fi
fi


 
