#!/bin/bash

dumpdir=$1
numfiles="$(($2-1))"
if [ $numfiles -lt 1 ]
then
    awk -f $SBUTIL/transpose_file.awk $dumpdir/dump0.csv > $dumpdir/dump0.dat
else
    i=0
    while [ $i -le $numfiles ]
    do
        awk -f $SBUTIL/transpose_file.awk $dumpdir/dump$i.csv > $dumpdir/dump$i.dat
	i="$(($i+1))"
    done
fi
