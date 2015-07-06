#!/bin/bash

for i in `cat ~/compact/reader/SC_SS0123_pass5_tot.list` ;
do
    #stager_qry -f ~/compact/reader/SC_SS0123_pass5_tot.list
    echo $i
    stager_get -M $i
    #if `grep -e STAGEIN` ; then
    #    stager_get -M $i
    #else
    #    if `grep -e Error`; then
    #        stager_get -M $i
    #    fi
    #fi

done
