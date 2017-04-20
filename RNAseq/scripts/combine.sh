#!/bin/bash
# combine.sh
# combines count data from all samples into one .txt file
# Seungsoo Kim

dir='../nobackup/counts'
paste <(cut -f 1 $dir/asy_r1.counts.txt) <(cut -f 2 $dir/asy_r1.counts.txt) <(cut -f 2 $dir/asy_r2.counts.txt) <(cut -f 2 $dir/asy_r3.counts.txt) <(cut -f 2 $dir/gal_r1.counts.txt) <(cut -f 2 $dir/gal_r2.counts.txt) <(cut -f 2 $dir/gal_r3.counts.txt) <(cut -f 2 $dir/sat_r1.counts.txt) <(cut -f 2 $dir/sat_r2.counts.txt) <(cut -f 2 $dir/sat_r3.counts.txt) > $dir/all.counts.txt
