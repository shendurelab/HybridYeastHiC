#!/bin/bash
# build_references.sh
# Seungsoo Kim

# load modules
module load bowtie2/2.2.3

# PARAMETERS
# directory for output
out="../nobackup"
if [ ! -d "$out" ]; then
  mkdir "$out"
fi
dir="../nobackup/annotations"
if [ ! -d "$dir" ]; then
  mkdir "$dir"
fi
# subdirectory for bowtie2 indices
bt2="../nobackup/bowtie2"
if [ ! -d "$bt2" ]; then
  mkdir "$bt2"
fi
# file with restriction enzymes - tab delimited file with following columns:
# 1) restriction enzyme name, 2) restriction site, 3) number of nt offset for cut
restriction_enzymes="../restr_enz.txt"
# file with references - tab-delimited file with following columns:
# 1) reference name, 2) reference fasta file path, 3) reference centromere annotation file path, 4) list of binsizes
refs="$1"

while read ref ref_file ref_cens binsizes
do
  echo "Building reference files for " $ref

  # make output folders
  for d in aligned assigned matrix
  do
    if [ ! -d "$out/$d/$ref" ]; then
      mkdir "$out/$d/$ref"
    fi
  done

  # build bowtie2 reference
  bowtie2-build ../$ref_file $bt2/$ref

  # digest genome with each restriction enzyme
  while read rname target offset
  do
    echo $rname
    ./digest_genome ../$ref_file $target $offset > $dir/$ref.$rname.bed
  done < ${restriction_enzymes}

  #create chrom_lengths file
  bowtie2-inspect -s $bt2/$ref | grep Sequence | grep -v '2micron' | grep -v 'MT' | cut -f2,3 > $dir/$ref.chrom_lengths 

  for binsize in $(echo $binsizes | tr "," "\n")
  do
    echo $binsize

    #make output folders
    if [ ! -d "$out/matrix/$ref/$binsize" ]; then
      mkdir "$out/matrix/$ref/$binsize"
    fi
  done
done < $refs
