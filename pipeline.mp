# A basic RNA Seq pipeline.
# Given a single sample, perform the following analysis.
# To analyze multiple samples, see the run script.
#
# author: Brian Schrader
# designer: Bhuvan Molparia
# since: 2016-03-31

[COMMANDS]
trimmomatic PE {1} {2} {o.fastq} {o.fastq} {o.fastq} {o.fastq} trim_options 
bowtie2 --very-sensitive -N 1 -p 8 -x HG_19 -q -1 {1.1-1} -2 {1.1-1} -S {o}
htseq-count {2.1} gene_list > {o.counts}

[PATHS]
trimmomatic java ../etc/Trimmomatic-0.36/trimmomatic-0.36.jar
bowtie2 ../etc/bowtie/bowtie2
htseq-count ../etc/htseq/htseq-count

HG_19 ../etc/hg_19
trim_options ILLUMINACLIP:Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10:2:true LEADING:3 TRAILING:3
gene_list ../etc/GENES

[FILES]
1. __read_1__
2. __read_2__
