#!/bin/sh
#
#$ -S /bin/bash
#$ -cwd
#$ -q jamaica
#$ -pe jamaica 8
#$ -N g_0-wtr-wtr
#

CUR=`pwd`
SRC=/pubdata/lhua/gms_"$JOB_ID"
mkdir -p $SRC
cp ../RunDimer.py $SRC
cp *.inp $SRC
cp *.inp.log $SRC
cd $SRC

python RunDimer.py

mv * $CUR && rm -rf $SRC
