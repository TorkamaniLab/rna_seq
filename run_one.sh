#! /bin/bash
# Runs the RNA Seq pipeline for a given sample.
# 
# Usage: ./run_one.sh <sample_id> <read_1> <read_2>
set -e
SAMPLE=$1
R1=$2
R2=$3

mkdir $SAMPLE;
cat pipeline.mp | 
    sed "s/__sample_id__/$SAMPLE/g" | 
    sed "s/__read_1_filename_/$R1/g" | 
    sed "s/__read_2_filename_/$R2/g" > $SAMPLE"/pipeline.mp"
cd $SAMPLE;
metapipe -j pbs -o pipeline.sh pipeline.mp
qsub pipeline.sh
