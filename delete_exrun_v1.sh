

#inputtag=k2pigd_kloe_pi0d_radcor
#inputtag=ee_split_2004
#inputtag=ee_split_2003
#inputtag=ee_pass5

inputtag=$1
tag=$2	
outversion=$3

lists=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/lists/$inputtag/${tag}_$outversion
outlogs=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/logs/$inputtag/${tag}_$outversion
output=/afs/cern.ch/user/r/rmarchev/compact/reader/tmp/output/$inputtag/${tag}_$outversion

rm -rf $lists $outlogs $output

 
