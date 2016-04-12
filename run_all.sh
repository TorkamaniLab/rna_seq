#! /bin/bash
# Given a list of samples, run them all.
# 
# author: Brian Schrader
# since: 2016-04-11
#
# Usage: ./run_all.sh <sample file>

SAMPLES=$1

cat $SAMPLES | while read SAMPLE FILE1 FILE2; do 
    bash run_one.sh $SAMPLE $FILE1 $FILE2;
done;

