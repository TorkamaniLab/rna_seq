# A basic RNA Seq pipeline.
# Given a single sample, perform the following analysis.
# To analyze multiple samples, see the run script.
#
# author: Brian Schrader, Bhuvan Molparia
# designer: Bhuvan Molparia
# since: 2017-09-21

[COMMANDS]
$TRIMMOMATIC PE {1} {2} {o.fastq} {o.fastq} {o.fastq} {o.fastq} $TRIM_OPTS 

#PBS -l npus=8
#PBS -l walltime=100:00:00
$BOWTIE2 --very-sensitive -N 1 -p 8 -x BOWTIE_INDEX -q -1 {1.1-1} -2 {1.1-3} -S {o}

$HTSEQ-COUNT -s no -t exon -i gene_id {2.1} $GTF_FILE > {o.counts}

[PATHS]
## Garibaldi paths 
$TRIMMOMATIC module load trimmomatic; java -jar `which trimmomatic-0.35.jar`
$BOWTIE2 module load bowtie2; bowtie2
$HTSEQ_COUNT module load python; htseq-count

## Local paths
#$TRIMMOMATIC /path/to/trimmomatic
#$BOWTIE2 /path/to/bowtie_2
#$HTSEQ_COUNT /path/to/htseq_count

$BOWITE_INDEX ../etc/bowtie_index
$GTF_FILE ../etc/gtf_file
$TRIM_OPTS ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3

[JOB_OPTIONS]
#PBS -l mem=12gb
set -e;

[FILES]
1. __read_1__
2. __read_2__
