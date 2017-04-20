#!/bin/bash
# process_reads.sh
# Seungsoo Kim

# load modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs
module load bowtie2/2.2.3 samtools/1.2

# PARAMETERS
# directory with sequencing data
dir="../data"
# directory for output
out="../nobackup"
# directory with references
refs="../references"
# directory with annotations
anns="../nobackup/annotations"
# subdirectory for bowtie2 indices
bt2="../nobackup/bowtie2"
# file with restriction enzymes - tab delimited file with following columns:
# 1) restriction enzyme name, 2) restriction site, 3) number of nt offset for cut
restriction_enzymes="../restr_enz.txt"

# optional parameters to override or avoid overriding existing files
override_map=true
override_asn=true

# sample name
samp=$1
# sample type (single or dual)
stype=$2
# reference name
ref=$3
# reference 1
ref1=$4
# reference 2
ref2=$5
# restriction enzyme name
renz=$6

echo $samp

# map reads separately
for i in 1 2
do
  echo "Mapping read" $i
  if [ $stype = "single" ]; then
    if [ ! -r $out/aligned/$ref/$samp.$i.bam ] || $override_map; then
      bowtie2 --very-sensitive -p $NSLOTS --reorder -x $bt2/$ref -U $out/junc_trim/$samp/${samp}_R$i.fastq.gz 2> $out/aligned/$ref/$samp.$i.bt2.out | samtools view -b - > $out/aligned/$ref/$samp.$i.bam
    fi
  else
    if [ ! -r $out/aligned/$ref1/$samp.$i.bam ] || $override_map; then
      bowtie2 --very-sensitive -p $NSLOTS --reorder -x $bt2/$ref1 -U $out/junc_trim/$samp/${samp}_R$i.fastq.gz 2> $out/aligned/$ref1/$samp.$i.bt2.out | samtools view -b - > $out/aligned/$ref1/$samp.$i.bam
    fi
    if [ ! -r $out/aligned/$ref2/$samp.$i.bam ] || $override_map; then
      bowtie2 --very-sensitive -p $NSLOTS --reorder -x $bt2/$ref2 -U $out/junc_trim/$samp/${samp}_R$i.fastq.gz 2> $out/aligned/$ref2/$samp.$i.bt2.out | samtools view -b - > $out/aligned/$ref2/$samp.$i.bam
    fi
  fi
done

echo "Assigning read pairs to restriction fragments"
rsite=$(awk -v r=$renz '{if ($1==r) print $2}' $restriction_enzymes)
if [ ! -r $out/assigned/$ref/$samp.out ] || $override_asn; then
  if [ $stype = "single" ]; then
    ./process_pairs $rsite $anns/$ref.$renz.bed 30 1000 <(samtools view $out/aligned/$ref/$samp.1.bam) <(samtools view $out/aligned/$ref/$samp.2.bam) $out/assigned/$ref/$samp.assigned $out/assigned/$ref/$samp.out
  else
    ./process_pairs_dual $rsite $anns/$ref.$renz.bed 30 1000 <(samtools view $out/aligned/$ref1/$samp.1.bam) <(samtools view $out/aligned/$ref2/$samp.1.bam) <(samtools view $out/aligned/$ref1/$samp.2.bam) <(samtools view $out/aligned/$ref2/$samp.2.bam) $out/assigned/$ref/$samp.assigned $out/assigned/$ref/$samp.out
  fi
  awk -v a=$ref1'_13' b=$ref2'_13' '($1==a && $7==b) || ($1==b && $7==a)' $out/assigned/$ref/$samp.chr13.assigned
fi
