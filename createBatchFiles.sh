#!/bin/bash

##############################################
#  Dis is a simple bash script which creates #
#  .sh files for LSF submition               #
##############################################

type="data"    #local folder where will be stores batch files
version="1"   #version of compact

#folder="data"  # run list files folder
folder=$type;
dir="/home/user/compact/reader/runlists/scmp36/period5/$folder"
mkdir -p "LSF/$version/$type"
mkdir -p "/gpfs/fs3/na62/user/user/LSF/$version/$type"


echo "Local folder $version/$type created"


ls "$dir" | while read -r file; 
do 
    bfile="LSF/$version/$type/$file.sh"
    tmpPath=""
    #echo "#!/bin/tcsh">> $bfile
    #echo "set tmpPath = $tmpPath" >>$bfile
    #echo "mkdir -p StmpPath">>$bfile
    echo 'set HOME=$PWD' >> $bfile
    echo "source /home/user/compact/reader/var.sh">>$bfile	
    echo "cd /home/mariov/compact/reader/">>$bfile
    echo "LSF=/jobdir/\${LSB_JOBID}">>$bfile	
    echo "/home/mariov/compact/reader/./compact  -l $dir/$file -rootfilename \$LSF/$version.$type.$file" >>$bfile
   # echo 'cd $LSF' >> $bfile
   # echo "xrdcp \$LSF/$file.root  $xroot$castordir/" >>$bfile;
    echo "cp \$LSF/$version.$type.$file.root /gpfs/fs3/na62/user/mariov/LSF/$version/$type/" >> $bfile
    echo "rm \$LSF/$version.$type.$file.root">>$bfile;
    chmod 755 $bfile
done

