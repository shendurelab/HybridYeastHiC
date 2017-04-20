#!/bin/bash
# pipeline.sh
# runs all scripts needed to analyze ChIP-seq data for one sample
# Seungsoo Kim

# load base modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs

# use cutadapt to trim adapters and low-quality ends of reads
module load python/2.7.3 cutadapt/1.8.3

out='../nobackup/cutadapt'
in='../data'
cutadapt -q 20 -m 28 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT -o $out/$1.1.fastq.gz -p $out/$1.2.fastq.gz $in/$2 $in/$3 > $out/$1.log

# use bowtie2 to align reads
module load htslib/1.3.2 bowtie2/2.2.3 samtools/1.3.1_htslib-1.3.2

in='../nobackup/cutadapt'
out='../nobackup/aligned'
bt2='../nobackup/bowtie2'

bowtie2 --very-sensitive -X 2000 -p $NSLOTS -x $bt2/sacCer3 -1 $in/$1.1.fastq.gz -2 $in/$1.2.fastq.gz 2> $out/$1.log | samtools view -b - > $out/$1.bam

# insert sizes
samtools view $out/$1.bam | awk '$5 >= 30 && $9 > 0 && $9 <= 2000 {OFS="\t"; print $9}' | head -n 10000 > $out/$1.insertsizes.txt

module load bedtools/2.24.0

in='nobackup/aligned'
out='nobackup/bedgraph'
genes='nobackup/references/sacCer3_genes.gff'
genome='nobackup/references/sacCer3.genome'

# convert BAM to BEDPE
bedtools bamtobed -i $in/$1.bam -bedpe > $in/$1.bedpe

# filter fragments <= 2000 bp with MAPQ >= 30 for both reads
awk '($1==$4) && ($6-$2 <= 2000) && ($8 >= 30) {OFS="\t"; print $1, $2, $6}' $in/$1.bedpe | sort -k1,1 -k2,2 -k3,3 > $in/$1.bed

# remove PCR duplicates
uniq $in/$1.bed > $in/$1.dedup.bed

# calculate genome coverage using bedtools
bedtools genomecov -i $in/$1.dedup.bed -g $genome -bg -trackline > $out/$1.bedgraph
