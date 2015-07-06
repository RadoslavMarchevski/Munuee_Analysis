
fnumber1=$1
fnumber2=$2

1=pass5_p
2=pass5_c_1st
3=pass5_c_2nd
4=pass5_c_3rd
5=pass5_c_4th
6=pass5_c_5th
7=pass5_c_6th

for i in `seq $fnumber1 $fnumber2` ;
do
    $i
    
     
    
    sh check_dead.sh $i
done
