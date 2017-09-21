#!/bin/bash

BOWTIE_INDEX=$1
GTF_FILE_PATH=$2

ln -T $BOWTIE_INDEX ./etc/bowtie_index
ln -T $GTF_FILE_PATH ./etc/gtf_file
