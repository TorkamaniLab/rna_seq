# A basic RNA Seq pipeline.
# Given a single sample, perform the following analysis.
# To analyze multiple samples, see the run script.
#
# author: Brian Schrader
# designer: Bhuvan Molparia
# since: 2016-03-31

[COMMANDS]
trimmomatic PE {1} {2} {o.fastq} {o.fastq} {o.fastq} {o.fastq} ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3
bowtie2 --very-sensitive -N 1 -p 8 -x HG_19 -q -1 {1.1-1} -2 {1.1-1} -S {o}
htseq-count {2.1} gene_list > {o.counts}
head {3.1} > summary.txt

[PATHS]
trimmomatic java etc/Trimmomatic-0.35/trimmomatic-0.35.jar
bowtie2 etc/bowtie/bowtie2
htseq-count etc/htseq/htseq-count

HG_19 etc/hg_19

[FILES]
1. __read_1__
2. __read_2__
