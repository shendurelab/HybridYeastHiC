#!/bin/bash
# run.sh
# main script, runs all scripts needed to analyze RNA-seq data
# Seungsoo Kim

# load base modules
. /etc/profile.d/modules.sh
module load modules modules-init modules-gs

# make output folders
for fold in 'cutadapt' 'aligned' 'counts' 'bedgraph' 'bowtie2'
do
  if [ ! -d '../nobackup/'$fold ]; then
    mkdir '../nobackup/'$fold
  fi
done

# make bowtie2 reference
if [ ! -f '../nobackup/bowtie2/sacCer3.bt2' ]; then
  module load bowtie2/2.2.3
  bowtie2-build ../nobackup/references/sacCer3.fa ../nobackup/bowtie2/sacCer3
fi


# run samples through pipeline
while read samp read1 read2
do
    qsub -N $samp -o ../nobackup/sge -e ../nobackup/sge -l mfree=2G -pe serial 4 -cwd ./pipeline.sh $samp $read1 $read2
done < samples.txt
