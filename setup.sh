#! /bin/bash
# Sets up the environment for the pipeline.
#
# Usage: ./setup.sh
set -e;

mkdir etc && cd etc;

# Download trimmomatic
wget -O trimmomatic.zip http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip;
unzip trimmomatic.zip
rm trimmomatic.zip


