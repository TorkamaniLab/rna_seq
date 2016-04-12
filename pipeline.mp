# A basic RNA Seq pipeline.
# Given a single sample, perform the following analysis.
# To analyze multiple samples, see the run script.
#
# author: Brian Schrader
# designer: Bhuvan Molparia
# since: 2016-03-31

[COMMANDS]
$TRIMMOMATIC PE {1} {2} {o.fastq} {o.fastq} {o.fastq} {o.fastq} $TRIM_OPTS 
$BOWTIE2 --very-sensitive -N 1 -p 8 -x HG_19 -q -1 {1.1-1} -2 {1.1-3} -S {o}
$HTSEQ-COUNT {2.1} $GTF_FILE > {o.counts}

[PATHS]
## Uncomment the versions you want.
## Garibaldi
$TRIMMOMATIC module load trimmomatic; trimmomatic 
$BOWTIE2 module load bowtie2; bowtie2
$HTSEQ_COUNT module load python; htseq-count

## Local
#$TRIMMOMATIC java ../etc/Trimmomatic-0.36/trimmomatic-0.36.jar
#$BOWTIE2 ../etc/bowtie/bowtie2
#$HTSEQ_COUNT ../etc/htseq/htseq-count

$HG_19 ../etc/hg_19
$TRIM_OPTS ILLUMINACLIP:Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3
$GTF_FILE ../etc/GENES

[JOB_OPTIONS]
#PBS -l mem=12gb
set -e;

[FILES]
1. __read_1__
2. __read_2__
