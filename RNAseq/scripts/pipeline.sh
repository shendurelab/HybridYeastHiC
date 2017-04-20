#!/bin/bash
# pipeline.sh
# runs all scripts needed to analyze RNA-seq data for one sample
# Seungsoo Kim

# load base modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs

# cutadapt - quality-trim and adapter-trim
module load python/2.7.3 cutadapt/1.8.3

out='../nobackup/cutadapt'
in='../data'

cutadapt -q 20 -m 28 -a AGATCGGAAGAGCACACGTC -A AGATCGGAAGAGCGTCGTGT -o $out/$1.1.fastq.gz -p $out/$1.2.fastq.gz $in/$2 $in/$3 > $out/$1.log

module load htslib/1.3.2 bowtie2/2.2.3 samtools/1.3.1_htslib-1.3.2

in='../nobackup/cutadapt'
out='../nobackup/aligned'
bt2='../nobackup/bowtie2'

# bowtie2 - align reads to sacCer3
bowtie2 --very-sensitive -p $NSLOTS -x $bt2/sacCer3 -1 $in/$1.1.fastq.gz -2 $in/$1.2.fastq.gz 2> $out/$1.log | samtools view -b - > $out/$1.bam

# sort BAM for genome coverage
samtools sort -o $out/$1.sorted.bam -T $out/$1 -@ $NSLOTS $out/$1.bam

# sample insert size distribution
samtools view $out/$1.bam | awk '$5 >= 30 && $9 > 0 && $9 <= 2000 {print $9}' | head -n 10000 > $out/$1.insertsizes.txt

module load bedtools/2.24.0 numpy/1.8.1 six/1.10.0 matplotlib/1.3.1 pysam/0.8.4 HTSeq/0.6.1p1

in='../nobackup/aligned'
out='../nobackup/counts'
genes='../nobackup/references/sacCer3_genes.gff'
genome='../nobackup/references/sacCer3.genome'

# count read pairs in each gene
python -m HTSeq.scripts.count -f bam -r name -s no -a 30 -t gene -i ID $in/$1.bam $genes > $out/$1.counts.txt

# use bedtools to generate genome coverage track
bedtools genomecov -ibam $in/$1.sorted.bam -g $genome -bg -trackline -pc > nobackup/$1.bedgraph
