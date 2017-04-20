#!/bin/bash
# create_Y12xDBVPG6044_mask.sh
# Seungsoo Kim

ref="Y12xDBVPG6044"
lib1="Y12_exponential_Sau3AI"
ref1="Y12"
lib2="DBVPG6044_exponential_Sau3AI"
ref2="DBVPG6044"

dir="../nobackup/aligned"

# load modules
module load samtools/1.2 bedtools/2.24.0

bedtools merge -i \
  <(paste <(cat <(samtools view $dir/$ref1/$lib1.1.bam) \
                <(samtools view $dir/$ref1/$lib1.2.bam) \
                <(samtools view $dir/$ref2/$lib2.1.bam) \
                <(samtools view $dir/$ref2/$lib2.2.bam)) \
          <(cat <(samtools view $dir/$ref2/$lib1.1.bam) \
                <(samtools view $dir/$ref2/$lib1.2.bam) \
                <(samtools view $dir/$ref1/$lib2.1.bam) \
                <(samtools view $dir/$ref1/$lib2.2.bam)) | \
    grep -w -v "4" | \
    cut -f2,3,4,5,6,12,14,15,16,17,18,21,22,23,24,25,31,33,34,35,36,37 | \
    awk '$4 >= 30 && $15 >= 30' | \
    awk '$6 != $17' | grep -v "XS:i:" | sed 's/AS:i://g' | sed 's/XM:i://g' | sed 's/XO:i://g' | grep -v 'XN' | sed 's/M//g' | \
    awk '($7 >= 1 && $7+$8 >= 2 && $18 == 0 && $19 == 0) {OFS="\t"; print $13, $14, $14+$16}' | \
    sort -k1,1 -k2,2n -k3,3n) \
  -c 1 -o count | \
  awk '$4 > 5' > ../masks/mismapped.bed
