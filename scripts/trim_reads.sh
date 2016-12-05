#!/bin/bash
# trim_reads.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load fastqc/0.11.3
module load python/2.7.3 cutadapt/1.8.3

# PARAMETERS
# directory with sequencing data
dir="../data"
# directory for output
out="../nobackup"
# file with restriction enzymes - tab delimited file with following columns:
# 1) restriction enzyme name, 2) restriction site, 3) number of nt offset for cut
restriction_enzymes="../restr_enz.txt"

# forward adapter sequence
fadapt='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
# reverse adapter seqeunce
radapt='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
# quality filter for trimming
qfilt=20
# minimum read length after trimming
minlen=20

run_fastqc=false
override_trim=false

# sample name
samp=$1
# restriction enzyme name
renz=$2
# read 1 file
r1=$3
# read 2 file
r2=$4

echo $samp

if $run_fastqc; then
  echo "Running fastqc"
  fastqc $r1
  fastqc $r2
fi

echo "Trimming adapter sequences"
# create output directory
if [ ! -d $out/adapt_trim/$samp ]; then
  mkdir $out/adapt_trim/$samp
fi
# run cutadapt to quality trim and trim adapter sequences
if [ ! -r $out/adapt_trim/$samp.adapt_trim.out ] || $override_trim; then
  cutadapt -q $qfilt -m $minlen -a $fadapt -A $radapt -o $out/adapt_trim/$samp/${samp}_R1.fastq.gz -p $out/adapt_trim/$samp/${samp}_R2.fastq.gz $dir/$r1 $dir/$r2 > $out/adapt_trim/$samp.adapt_trim.out
fi

echo "Trimming to ligation junction"
# determine ligation junction sequence
junc=$(awk -v r=$renz '{if ($1==r) print substr($2,1,length($2)-$3) substr($2,$3+1,length($2)-$3)}' $restriction_enzymes)
# create output directory
if [ ! -d $out/junc_trim/$samp ]; then
  mkdir $out/junc_trim/$samp
fi
# run cutadapt to trim reads to ligation junction
if [ ! -r $out/junc_trim/$samp.junc_trim.out ] || $override_trim; then
  cutadapt -m $minlen -a $junc -A $junc -o $out/junc_trim/$samp/${samp}_R1.fastq.gz -p $out/junc_trim/$samp/${samp}_R2.fastq.gz $out/adapt_trim/$samp/${samp}_R1.fastq.gz $out/adapt_trim/$samp/${samp}_R2.fastq.gz > $out/junc_trim/$samp.junc_trim.out
fi
