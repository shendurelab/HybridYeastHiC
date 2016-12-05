#!/bin/bash
# make_matrix.sh
# shell script for make_matrix.cpp
# Seungsoo Kim

dir="../references"
anns="../nobackup/annotations"
out="../nobackup"
masks="../masks"

samp=$1
ref=$2
bsize=$3
mask=$4
minavg=$5

./make_matrix $anns/$ref.chrom_lengths $bsize $masks/$mask.bed <(grep -v invalid $out/assigned/$ref/$samp.assigned) $minavg $out/matrix/$ref/$bsize/$samp.$mask.rowsums $out/matrix/$ref/$bsize/$samp.$mask.matrix
