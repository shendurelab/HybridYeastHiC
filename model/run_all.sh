#!/bin/bash
# run_all.sh
# 
# Shell script to run multiple iterations of structure simulations

# compile C++ script to process PDB output files
g++ pdb2interactions.cpp -o pdb2interactions -O3

# make SGE log folders
if [ ! -d nobackup ]; then
  mkdir nobackup
fi
if [ ! -d nobackup/sge_out ]; then
  mkdir nobackup/sge_out
fi
if [ ! -d nobackup/sge_error ]; then
  mkdir nobackup/sge_error
fi

chrlens='../nobackup/annotations/ScSu.chrom_lengths'
cens='../references/ScSu_2micron_revised_cen.gff'

# run 20,000 iterations of model for S. cerevisiae x S. uvarum hybrid for full diploid size (1.25x each dimension compared to haploid model)
if [ ! -d ScSu_v1 ]; then
  mkdir ScSu_v1
fi
qsub -t 1-20000 -q ravana.q -N ScSu_v1 -o nobackup/sge_out/ -e nobackup/sge_error/ -l mfree=1G -cwd -b y -tc 80 ./run.sh ScSu_v1 $chrlens $cens 1.25

# run 20,000 iterations of model for S. cerevisiae x S. uvarum hybrid for 80% size
if [ ! -d ScSu_v2 ]; then
  mkdir ScSu_v2
fi
qsub -t 1-20000 -q ravana.q -N ScSu_v2 -o nobackup/sge_out/ -e nobackup/sge_error/ -l mfree=1G -cwd -b y -tc 80 ./run.sh ScSu_v2 $chrlens $cens 1.0

# run 20,000 iterations of model for S. cerevisiae x S. uvarum hybrid for 64% size
if [ ! -d ScSu_v3 ]; then
  mkdir ScSu_v3
fi
qsub -t 1-20000 -q ravana.q -N ScSu_v3 -o nobackup/sge_out/ -e nobackup/sge_error/ -l mfree=1G -cwd -b y -tc 80 ./run.sh ScSu_v3 $chrlens $cens 0.8

