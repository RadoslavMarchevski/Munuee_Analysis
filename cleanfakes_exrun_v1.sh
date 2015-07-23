

#inputtag=k2pigd_kloe_pi0d_radcor
#inputtag=ee_split_2004
#inputtag=ee_split_2003
#inputtag=ee_pass5

inputtag=$1
tag=$2	
outversion=$3
size=$4


output=/home/rmarchev/compact/reader/tmp/output/$inputtag/${tag}_$outversion


for i in $output/[1-9]*/* ; do  
	FILESIZE=$(stat -c%s "$i");  
	if [ $FILESIZE -le $size ]; then 
		echo $i	"	" $FILESIZE; 
		#rm $i; 
	fi; 
done


 
