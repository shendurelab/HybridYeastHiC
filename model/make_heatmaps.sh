#!/bin/bash
# make_heatmaps.sh
# 
# Processes interaction counts to make Hi-C heatmaps (matrices)
# Seungsoo Kim

# creates list corresponding each bead number to 32 kb bins, using PDB file
for ref in ScSb_v1 ScSb_v2 ScSb_v3
do
  awk 'BEGIN{OFS="\t"; chr=""; b=0; cum=-1}{if(chr != $10) cum+=int((b-1)/10)+1; print cum+int(($6-1)/10), $10, $3, $6; b=$6; chr=$10}' $ref/1.pdb > $ref.32000.bins
done

# compile C++ script to aggregate interaction counts into bins, to make heatmaps
g++ interactions2heatmap.cpp -o interactions2heatmap -O3

# run C++ script through simulation data
for ref in 'ScSb_v1' 'ScSb_v2' 'ScSb_v3'
do
  for c in 45
  do
    ./interactions2heatmap $ref.32000.bins $c 1.0 <(cat $ref/{1..20000}.c$c.interactions.txt) $ref.c$c.32000.rowsums $ref.c$c.32000.matrix
  done
done
