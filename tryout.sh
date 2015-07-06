#!/bin/bash

HOME=$PWD
echo "Hello World"

LSF=${LSB_CWD}/${LSB_JOBID}
echo "$LSF"
cd /afs/cern.ch/user/r/rmarchev/compact/reader/

make

#/afs/cern.ch/user/r/rmarchev/compact/reader/./compact -l tail.list  -rootfilename ${LSF}/tryout -cheat 

#cp -f ${LSF}/tryout.root cd /afs/cern.ch/user/r/rmarchev/compact/reader/test/tryout.root
